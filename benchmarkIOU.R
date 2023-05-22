library(patchwork)
library(ggpubr)
library(extrafont)

library(microbenchmark)
library(bit)
library(tidyverse)
library(magrittr)
library(conflicted)
c("symdiff", "xor") %>% 
  sapply(function(x) conflict_prefer(x, "bit", quiet = T)) %>% 
  invisible()
conflict_prefer("filter", "tidygraph")
c("group_rows", "filter", "select", "lag") %>% 
  sapply(function(x) conflict_prefer(x, "dplyr", quiet = T)) %>% 
  invisible()
c("set_names", "set_colnames", "set_rownames") %>% 
  sapply(function(x) conflict_prefer(x, "magrittr", quiet = T)) %>% 
  invisible()
c("expand", "pack", "unpack") %>% 
  sapply(function(x) conflict_prefer(x, "tidyr", quiet = T)) %>% 
  invisible()

source("IOUHelpers.R")
liou <- function(l1, l2) {
  sum(l1 & l2) / sum(l1 | l2)
}
bit_vec_s3iou <- function(v1, v2) {
  sum(v1 & v2) / sum(v1 | v2)
}

downscale_person <- function(p, d) {
  p$start <- ceiling(p$start / d)
  p$end   <- floor(p$end / d)
  p$length <- floor(p$length / d)
  if (is.na(p$length)) stop("Length of person after rounding is NA/NaN")
  p
}

# Before running this run the Hyperparameter cell in the Rmd with:
# - d = 1
# - segment_type = "all"
# - chromosome = "all"

set.seed(123)

tprs_ind <- sample(nrow(arch_persons), 2)
tprs_1 <- arch_persons$prs[[tprs_ind[1]]]
tprs_2 <- arch_persons$prs[[tprs_ind[2]]]


benchmark_results <- tibble(
  downscale = 10^seq(0, 4, .1)
) %>% 
  mutate(
    res = map(downscale, function(d) {
      p1 <- downscale_person(tprs_1, d)
      p2 <- downscale_person(tprs_2, d)
      l1 <- do.call(encode_logical_vector, p1)
      l2 <- do.call(encode_logical_vector, p2)
      b1 <- do.call(encode_binary_vector, p1)
      b2 <- do.call(encode_binary_vector, p2)
      
      microbenchmark(
        implicit = implicit_iou(p1, p2),
        bit      = bit_vec_s3iou(b1, b2),
        bitCustom= bit_vec_iou(b1, b2),
        logical  = liou(l1, l2),
        times = 10
      ) %>% 
        as_tibble
    },
    .progress = T)
  ) %>% 
  unnest(res)


benchmark_plot <- benchmark_results %>% 
  mutate(
    expr = c(
       "implicit" = "implicit",
       "bit" = "bit (S3)",
       "bitCustom" = "bit (custom)",
       "logical" = "base"
    )[expr]
  ) %>% 
  ggplot(aes(downscale, time / 10^6, color = expr)) +
  stat_summary(geom = "line",
               fun = median,
               linewidth = .75,
               key_glyph = draw_key_point) +
  scale_y_log10(labels = prettify_scientific) +
  scale_x_log10() +
  scale_color_brewer(palette = "Set2") +
  coord_cartesian(expand = F) +
  guides(color = guide_legend(override.aes = list(size = 5,
                                                  shape = 16))) +
  labs(x = "Downscaling factor", y = "Time (ms)", color = "Method") +
  theme(
    legend.position = c(.25, .25),
    legend.box.background = element_rect(color = "black",
                                         linetype = "solid",
                                         linewidth = 1),
    plot.margin = margin(0,20,0,2)
  )

ggsave("benchmark.pdf", benchmark_plot,
       device = cairo_pdf, antialias = "subpixel",
       width = 4, height = 3.7, scale = 1.5, dpi = 200)
