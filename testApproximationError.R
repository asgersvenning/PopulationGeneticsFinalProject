library(patchwork)
library(ggpubr)
library(extrafont)
library(ggtext)

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
source("miscHelpers.R")

create_breaks <- function(len, lambda, freq) {
  brs <- list(start = c(), end = c(), length = len)
  sec_start <- 1
  last_also <- F
  while (TRUE) {
    sec_len   <- rpois(1, lambda = lambda)
    sec_end   <- sec_start + sec_len
    if (sec_end >= len) break
    if (rbinom(1, 1, freq) == 1) {
      if (length(brs$end) != 0 && (sec_start - 1) == brs$end[length(brs$end)]) {
        brs$end[length(brs$end)] <- sec_end
      }
      else {
        brs$end   <- c(brs$end, sec_end)
        brs$start <- c(brs$start, sec_start)
      }
    }
    sec_start <- sec_end + 1
  }
  
  brs
}

bit_vec_iou_split <- function(p1, p2, split) {
  split_person <- function(p, s) {
    list(
      list(
        start = p$start[p$start <= s],
        end   = p$end[p$end <= s],
        length = s
      ),
      list(
        start = p$start[p$start > s] - s,
        end   = p$end[p$end > s] - s,
        length = p$length - s
      )
    )
  }
  
  p1 <- split_person(p1, split)
  p2 <- split_person(p2, split)
  print(str(p1))
  print(str(p2))
  
  p1 <- lapply(p1, function(x) do.call(encode_binary_vector, x))
  p2 <- lapply(p2, function(x) do.call(encode_binary_vector, x))
  int1 <- sum(p1[[1]] & p2[[1]])
  int2 <- sum(p1[[2]] & p2[[2]])
  uni1 <- sum(p1[[1]] | p2[[1]])
  uni2 <- sum(p1[[2]] | p2[[2]])
  
  (int1 + int2) / (uni1 + uni2)
} 

# seed = 123, is used for the report (seed = 135 for two Papuans)
set.seed(123)
# tprs_1 <- create_breaks(10^7, 1000, .1)
# tprs_2 <- create_breaks(10^7, 10, .1)
tprs_ind <- sample(nrow(arch_persons), 2)
tprs_1 <- arch_persons$prs[[tprs_ind[1]]]
tprs_2 <- arch_persons$prs[[tprs_ind[2]]]
arch_persons %>% 
  slice(tprs_ind)

split_ind <- 2^31 - 1
stopifnot(sum(map2_lgl(tprs_1$start, tprs_1$end, function(s,e) between(split_ind, s, e))) == 0)
stopifnot(sum(map2_lgl(tprs_2$start, tprs_2$end, function(s,e) between(split_ind, s, e))) == 0)

true_iou <- bit_vec_iou_split(tprs_1, tprs_2, split_ind)
print(true_iou)

approxError <- tibble(
  downscale = 10^seq(3, 6, .05)
) %>% 
  mutate(
    iou = map(downscale, function(x) {
      testTypes <- tibble(type = c(rep("random", 5), "round", "floor", "ceiling", "expand", "contract")) 
      
      testTypes %>% 
        mutate(
          iou = map_dbl(type, ~bit_vec_iou(create_bit_vec(tprs_1, x, .x),  create_bit_vec(tprs_2, x, .x)))
        )
    }, .progress = T)
  ) %>% 
  unnest(iou) %>% 
  mutate(res = (iou - true_iou)/true_iou)

approxErrorPlt <- approxError %>% 
  group_by(downscale, type) %>% 
  summarize(
    res = mean(res),
    iou = mean(iou),
    .groups = "drop"
  ) %>% 
  mutate(
    type = factor(snakecase::to_sentence_case(type), c(
      "Ceiling",
      "Floor",
      "Round",
      "Random",
      "Expand",
      "Contract"
    ))
  ) %>% 
  filter(downscale <= 100000 & downscale != 1) %>%
  ggplot(aes(downscale, res, color = type)) +
  geom_line(linewidth = .75,
            key_glyph = draw_key_point) +
  geom_hline(yintercept = 0,#approxError$iou[approxError$downscale == 1][1],
             color = "black", linetype = "dashed", linewidth = .5) +
  scale_x_log10(labels = prettify_scientific, expand = expansion(),
                n.breaks = 10) +
  scale_y_continuous(labels = scales::label_percent(1),
                     trans = abslog(mult = 20),
                     breaks = c(-2, -1, -.5, -.2, -.1, -.05, -.01, 
                                # 0, 
                                .01, .05, .1, .2, .5, 1, 2, 5, 10),
                     limits = c(-.1, 2)) +
  scale_color_brewer(palette = "Dark2") +
  coord_cartesian(expand = F) +
  guides(color = guide_legend(override.aes = list(shape = 16,
                                                  size = 5,
                                                  alpha = 1),
                              ncol = 3)) +
  ggpubr::theme_pubr(legend = "right", base_family = "CMU Serif") +
  labs(color = "Segment approximation\nmethod",
       # caption = "<b><span style='font-size: 12pt;'>Segments are approximated by dividing the start & end indices by a constant, shown on the x-axis, followed by one of the following rounding methods:</span></b><br>- <b>Ceiling:</b> both start and end indices are rounded up.<br>- <b>Floor:</b> both start and end indices are rounded down.<br>- <b>Round</b>: both start and end indices are rounded as normal.<br>- <b>Random:</b> both start and end indices are \"random-rounded\", i.e. rounded up or down by a coin-flip weighted by their decimal part.<br>- <b>Expand:</b> the start indices are rounded down, while the end indices are rounded up.<br>- <b>Contract:</b> the start indices are rounded up, while the end indices are rounded down.",
       y = "<b>Relative</b> d<sub>IoU </sub> <b>error</b>", x = "Downscaling factor") +
  theme(
    axis.title.y = element_textbox(orientation = "left-rotated",
                                   size = 18),
    axis.title.x = element_text(face = "bold",
                              size = 18),
    axis.text = element_text(size = 13),
    panel.grid.major.y = element_line(colour = "gray80", linetype = "dashed", linewidth = .25),
    legend.position = c(0.325, .75),
    legend.title = element_text(face = "bold",
                                size = 18,
                                hjust = .5),
    legend.text = element_text(size = 14),
    legend.box.background = element_rect(colour = "black",
                                         linewidth = 1),
    plot.caption = element_textbox_simple(hjust = 0,
                                          box.colour = "black",
                                          width = unit(40, "lines"),
                                          linetype = "solid",
                                          linewidth = .5,
                                          fill = "gray85",
                                          margin = margin(1,0,0,0,"lines"),
                                          padding = margin(.5,.5,.5,.5,"lines")),
    plot.margin = margin(0.5, 1, 0.1, 0.1, "lines")
  )


ggsave("approxError_zoom.pdf", approxErrorPlt,
       device = cairo_pdf,
       width = 4, height = 4, scale = 2, dpi = 400)
