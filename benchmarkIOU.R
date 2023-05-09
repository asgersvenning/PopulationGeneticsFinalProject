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
  sum(l1 & l2) / sum(l1 | l2)
}

t1 <- rbinom(10^8,1,.5) %>% as.bit()
t2 <- !t1
l1 <- as.logical(t1)
l2 <- as.logical(t2)
t3 <- rbinom(10^8,1,.5) %>% as.bit()
t4 <- rbinom(10^8,1,.5) %>% as.bit()
l3 <- as.logical(t3)
l4 <- as.logical(t4)

bench <- microbenchmark(
  bit_vec_iou(t1, t2),
  liou(l1, l2),
  bit_vec_iou(t3, t4),
  liou(l3, l4),
  times = 100L
)

bench_plt <- bench %>% 
  as_tibble()  %>% 
  mutate(type = ifelse(str_detect(expr, "liou"), "logical", "bit")) %>% 
  arrange(type, desc(expr)) %>% 
  mutate(
    expr = factor(expr, unique(expr))
  ) %>% 
  ggplot(aes(time/10^6, expr, fill = type)) +
  geom_violin(scale = "width",
              color = "black",
              key_glyph = draw_key_point) +
  scale_x_log10(limits = c(50, 2000),
                expand = expansion()) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Milliseconds",
       y = NULL,
       fill = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21,
                                                 size = 5,
                                                 color = "transparent")))


size_plt <- tibble(
  type = c("bit", "logical"),
  size = c(object.size(t1), object.size(l1))
) %>% 
  ggplot(aes(size/10^6, type, fill = type)) +
  geom_col(color = "black",
           key_glyph = draw_key_point) +
  scale_x_log10(limits = c(1, 1000),
                expand = expansion()) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Mb", y = NULL, title = "", fill = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21,
                                                 size = 5,
                                                 color = "transparent")))


list(bench_plt, size_plt) %>% 
  wrap_plots(ncol = 1, guides = "collect") &
  # plot_annotation(
  #   title = "array size = 10^8", 
  #   theme = theme(plot.title = element_text(size = 24, 
  #                                           hjust = .5, 
  #                                           family = "CMU Serif"))
  # ) &
  theme_pubr(legend = "right") +
  theme(title = element_text(face = "bold",
                             size = 14)) 



ts1_4 <- rbinom(10^4, 1, .5) %>% as.bit()
ts2_4 <- rbinom(10^4, 1, .5) %>% as.bit()
ls1_4 <- as.logical(ts1_4)
ls2_4 <- as.logical(ts2_4)
ts3_5 <- rbinom(10^5, 1, .5) %>% as.bit()
ts4_5 <- rbinom(10^5, 1, .5) %>% as.bit()
ls3_5 <- as.logical(ts3_5)
ls4_5 <- as.logical(ts4_5)
ts5_6 <- rbinom(10^6, 1, .5) %>% as.bit()
ts6_6 <- rbinom(10^6, 1, .5) %>% as.bit()
ls5_6 <- as.logical(ts5_6)
ls6_6 <- as.logical(ts6_6)
ts7_7 <- rbinom(10^7, 1, .5) %>% as.bit()
ts8_7 <- rbinom(10^7, 1, .5) %>% as.bit()
ls7_7 <- as.logical(ts7_7)
ls8_7 <- as.logical(ts8_7)
ts9_8 <- rbinom(10^8, 1, .5) %>% as.bit()
ts10_8 <- rbinom(10^8, 1, .5) %>% as.bit()
ls9_8 <- as.logical(ts9_8)
ls10_8 <- as.logical(ts10_8)
                  
bench_2 <- microbenchmark(
  bit_vec_iou(ts1_4, ts2_4),
  bit_vec_iou(ts3_5, ts4_5),
  bit_vec_iou(ts5_6, ts6_6),
  bit_vec_iou(ts7_7, ts8_7),
  bit_vec_iou(ts9_8, ts10_8),
  liou(ls1_4, ls2_4),
  liou(ls3_5, ls4_5),
  liou(ls5_6, ls6_6),
  liou(ls7_7, ls8_7),
  liou(ls9_8, ls10_8),
  times = 100L
)


bench_2 %>% 
  as_tibble() %>% 
  mutate(type = str_detect(expr, "^liou"),
         type = ifelse(type, "logical", "bit"),
         size = str_extract(expr, "(?<=_)\\d"),
         size = 10^as.integer(size)) %>% 
  ggplot(aes(size, time/10^9, color = type, fill = type)) +
  stat_summary(geom = "ribbon", 
               color = "transparent", 
               fun.data = median_hilow,
               alpha = .25) +
  stat_summary(geom = "path", linewidth = .75) +
  scale_x_log10(labels = Vectorize(function(x) if (is.na(x) || !str_detect(x, "e\\+")) return(x) else latex2exp::TeX(paste0("$", str_replace(x, "e\\+", "0^{"),"}$")))) +
  scale_y_log10(labels = Vectorize(function(x) if (is.na(x) || !str_detect(x, "e\\-")) return(x) else latex2exp::TeX(paste0("$", str_replace(x, "e\\-", "0^{-"),"}$")))) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1", guide = "none") +
  theme_pubr(base_family = "CMU Serif", legend = "right") +
  labs(x = "Size of vectors in bits", y = "Time in seconds") +
  coord_cartesian(expand = F)




ls_comp <- tibble(
  len = seq(0, 6, .01)
) %>% 
  mutate(
    len = floor(10^len)
  ) %>% 
  distinct(len) %>% 
  slice_sample(n = nrow(.), replace = F) %>% 
  mutate(
    lgl = map(len, function(x) {
      l1 <- rbinom(x, 1, .5) == 1
      l2 <- rbinom(x, 1, .5) == 1
      tibble(
        size = object.size(l1),
        time = microbenchmark(liou(l1, l2), times = 10L)$time %>% 
          mean %>% 
          divide_by(10^9)
      )
      }),
    bit = map(len, function(x) {
      b1 <- as.bit(rbinom(x, 1, .5) == 1)
      b2 <- as.bit(rbinom(x, 1, .5) == 1)
      tibble(
        size = object.size(b1),
        time = microbenchmark(bit_vec_iou(b1, b2), times = 10L)$time %>% 
          mean %>% 
          divide_by(10^9)
      )
      })
  )


prettify_scientific <- Vectorize(function(x) {
  if (is.na(x)) return(NA)
  if (!is.character(x)) x <- as.character(x)
  if (!str_detect(x, "e\\+|e\\-")) return(x)
  if (str_detect(x, "e\\+")) {
    latex2exp::TeX(paste0("$", str_replace(x, "1*e\\+", " \\\\cdot 10^{"),"}$"))
  }
  else if (str_detect(x, "e\\-")) {
    latex2exp::TeX(paste0("$", str_replace(x, "1*e\\-", " \\\\cdot 10^{-"),"}$"))
  }
  else {
    stop("Invalid format!")
  }
})

ls_comp %>% 
  pivot_longer(!len, names_to = "type", values_to = "res") %>%
  unnest(res) %>% 
  mutate(size = as.numeric(size)) %>% 
  pivot_longer(!c(len, type), names_to = "var") %>% 
  ggplot(aes(len,value,color=type)) +
  geom_line(linewidth = .75) +
  scale_x_continuous(labels = prettify_scientific) +
  scale_y_continuous(labels = prettify_scientific) +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~var, scales = "free_y") +
  labs(y = NULL, x = "Length of vector/array") +
  coord_cartesian(expand = F) +
  ggpubr::theme_pubr() +
  theme(title = element_text(face = "bold"),
        plot.margin = margin(0,1,.25,0.1, "cm"))
