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

tprs_1 <- create_breaks(10^6, 1000, .1)
tprs_2 <- create_breaks(10^6, 10, .1)
# tprs_1 <- arch_persons$prs[[1]]
# tprs_2 <- arch_persons$prs[[2]]

approxError <- tibble(
  downscale = c(1, 10^seq(1, 7, .05))
) %>% 
  mutate(
    iou = map(downscale, function(x) {
      testTypes <- tibble(type = c(rep("random", 5), "round", "floor", "ceiling", "expand", "contract")) 
      
      testTypes %>% 
        mutate(
          iou = map_dbl(type, ~bit_vec_iou(create_bit_vec(tprs_1, x, .x),  create_bit_vec(tprs_2, x, .x)))
        )
    })
  ) %>% 
  unnest(iou) %>% 
  mutate(res = (iou - iou[downscale == 1][1])/iou[downscale == 1][1])

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
  # filter(downscale <= 100000) %>%
  ggplot(aes(downscale, res, color = type)) +
  geom_line(linewidth = .75,
            key_glyph = draw_key_point) +
  geom_hline(yintercept = 0,#approxError$iou[approxError$downscale == 1][1],
             color = "black", linetype = "dashed", linewidth = .5) +
  scale_x_log10(labels = prettify_scientific, expand = expansion(),
                n.breaks = 10) +
  scale_y_continuous(labels = scales::label_percent()) +
  scale_color_brewer(palette = "Dark2") +
  guides(color = guide_legend(override.aes = list(shape = 16,
                                                  size = 5,
                                                  alpha = 1))) +
  ggpubr::theme_pubr(legend = "right", base_family = "CMU Serif") +
  labs(color = "Segment\napproximation",
       caption = "<b><span style='font-size: 12pt;'>Segments are approximated by dividing the start & end indices by a constant, shown on the x-axis, followed by one of the following rounding methods:</span></b><br>- <b>Ceiling:</b> both start and end indices are rounded up.<br>- <b>Floor:</b> both start and end indices are rounded down.<br>- <b>Round</b>: both start and end indices are rounded as normal.<br>- <b>Random:</b> both start and end indices are \"random-rounded\", i.e. rounded up or down by a coin-flip weighted by their decimal part.<br>- <b>Expand:</b> the start indices are rounded down, while the end indices are rounded up.<br>- <b>Contract:</b> the start indices are rounded up, while the end indices are rounded down.",
       y = "Relative IoU error", x = "Downscaling factor") +
  theme(
    axis.title = element_text(face = "bold",
                              size = 18),
    axis.text = element_text(size = 13),
    legend.position = c(0.25, .75),
    legend.title = element_text(face = "bold",
                                size = 18,
                                hjust = .5),
    legend.text = element_text(size = 14),
    legend.box.background = element_rect(colour = "black",
                                         linewidth = .5),
    plot.caption = element_textbox_simple(hjust = 0,
                                          box.colour = "black",
                                          width = unit(40, "lines"),
                                          linetype = "solid",
                                          linewidth = .5,
                                          fill = "gray85",
                                          margin = margin(1,0,0,0,"lines"),
                                          padding = margin(.5,.5,.5,.5,"lines")),
    plot.margin = margin(0.1, 1, 0.1, 0.1, "lines")
  )


ggsave("approxError_zoom.pdf", approxErrorPlt,
       device = cairo_pdf,
       width = 8, height = 6, scale = 2, dpi = 400)
