# Data wrangling
library(magrittr)
library(tidyverse) %>% 
  suppressPackageStartupMessages()

# Plotting
library(ggforce)
library(patchwork)
library(ggnewscale)
library(extrafont)
library(ggpubr)
library(kableExtra)

# Graph utils
library(tidygraph)
library(ggraph)

# Projection utils
library(uwot)
library(ape)

# Bit array library for super speed and micro-memory usage
library(bit)

# Progress bar
library(progress)

# Resolve some name space conflicts
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

# Helper functions for storing, manipulating and converting between efficient interval-based 
# representations of archaic segments and binary bit arrays, as well as ultra efficient calculation
# of bit array intersection over union (IoU)
source("IOUHelpers.R")

arch <- read_delim("ArchaicSegments.txt", delim = "\t") %>% 
  # filter(chrom == "3") %>%
  filter(length > 0) %>% 
  mutate(chrom = factor(chrom, unique(chrom)[order(as.numeric(unique(chrom)))]))

all_chrom_range <- arch %>%
  group_by(chrom) %>% 
  summarize(
    start = min(start),
    end = max(end),
    .groups = "drop"
  ) %>% 
  arrange(chrom) %>% 
  mutate(
    length = end - start,
    offset = lag(cumsum(length), default = 0) #- first(length)
  )

all_chrom_summary <- all_chrom_range %>% 
  summarize(across(!chrom, ~list(set_names(.x, chrom)))) %>% 
  as.list %>% 
  lapply(first)

all_chrom_range %>% 
  mutate(across(c(start, end), ~.x - start + 1 + offset)) %>% 
  mutate(chrom = as.integer(chrom)) %>%
  mutate(nchrom = lead(chrom),
         nstart = lead(start)) %>% 
  ggplot() +
  geom_segment(aes(chrom, start, xend = chrom, yend = end)) +
  geom_segment(aes(chrom, end, xend = nchrom, yend = nstart)) +
  scale_x_continuous(labels = levels(all_chrom_range$chrom), breaks = 1:length(levels(all_chrom_range$chrom)))


arch_approx <- arch %>% 
  # group_by(chrom) %>% 
  # mutate(across(c(start, end), ~.x - all_chrom_summary$start[first(chrom)] + all_chrom_summary$offset[first(chrom)] + 1)) %>%
  mutate(across(c(start, end), ~.x/10000),
         start = ceiling(start),
         end = floor(end)) %>%
  mutate(length = end - start) %>% 
  filter(length > 0)

x_chrom_end_range <- arch_approx %>%
  pull(end) %>% 
  range

x_chrom_start_range <- arch_approx %>% 
  pull(start) %>% 
  range

x_chrom_min <- x_chrom_start_range[1]
x_chrom_max <- x_chrom_end_range[2] - x_chrom_min + 1

arch_approx <- arch_approx %>%
  mutate(across(c(start, end), ~.x - x_chrom_min + 1))

arch_persons <- arch_approx %>% 
  mutate(country = str_extract(country, "\\w+")) %>% 
  group_by(name, pop, country, region) %>% 
  summarize(
    prs = list(person(start, end, x_chrom_max)),
    .groups = "drop"
  ) %>% 
  mutate(prs = prs %>% 
           set_names(name))

all_dmat <- arch_persons %>%
  # slice_sample(n = 5) %>%
  person_pairwise_dist_df("prs", "name")

### Create plot
skip_labels <- function(n) function(x) replace(x, which(!(1:length(x) %in% floor(seq(1, length(x), length.out = n)))), "")

prs_to_country <- all_dmat %>% 
  distinct(country_1, prs_1) %>% 
  {set_names(.$country_1, .$prs_1)}

label_prs_to_country <- function(prs) {
  country <- prs_to_country[prs] %>% 
    unname
  
  country_boundary_ind <- c(1, which(country[1:(length(country) - 1)] != country[2:length(country)]))
  
  country_center_ind <- country_boundary_ind + 
    (c(country_boundary_ind, length(country)) %>% 
       diff %>% 
       divide_by(2) %>% 
       ceiling)
  
  country[-country_center_ind] <- "" # str_extract(country[-country_center_ind], ".")
  # cat(paste0(country, collapse = "\n"))
  country
}

all_dmat_plt <- all_dmat %>% 
  # filter(country_1 == "Papua" & country_2 == "Papua") %>% 
  arrange(region_1, country_1, pop_1, prs_1) %>%  
  mutate(
    country_1 = factor(country_1, unique(country_1)) %>% 
      fct_shuffle,
    country_2 = factor(country_2, levels(country_1)),
    prs_1 = factor(prs_1, unique(prs_1)),
    prs_2 = factor(prs_2, levels(prs_1))
  ) %>% 
  mutate(iou = ifelse(prs_1 == prs_2, NA, iou)) %>% 
  ggplot(aes(prs_1, prs_2, fill = iou, color = iou)) +
  geom_tile() +
  scale_fill_viridis_c(#trans = "log10",
    option = "A",
    direction = -1,
    na.value = "gray75",
    labels = scales::label_percent(),
    name = "1 - IoU") +
  scale_colour_viridis_c(#trans = "log10",
    option = "A",
    direction = -1,
    na.value = "gray75",
    labels = scales::label_percent()) +
  new_scale_fill() +
  geom_tile(aes(x = -4.5, y = ifelse(as.integer(prs_1) == 1, prs_2, NA), fill = as.integer(country_2), color = after_scale(fill)),
            width = 10,
            show.legend = F) +
  geom_tile(aes(y = -4.5, x = ifelse(as.integer(prs_2) == 1, prs_1, NA), fill = as.integer(country_1), color = after_scale(fill)),
            height = 10,
            show.legend = F) +
  scale_fill_distiller(palette = "Greys") +
  scale_x_discrete(labels = label_prs_to_country, expand = expansion(0, c(10, 0))) +
  scale_y_discrete(labels = label_prs_to_country, expand = expansion(0, c(10, 0))) +
  coord_equal() +
  guides(colour = guide_none()) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(hjust = 1,
                                   angle = 90,
                                   vjust = .1),
        axis.text.y = element_text(vjust = .1),
        axis.text = element_text(size = 6),
        axis.title = element_blank(),
        legend.title = element_text(hjust = .5,
                                    size = 20),
        legend.text = element_text(size = 14),
        legend.key.width = unit(2, "lines"),
        legend.key.height = unit(4, "lines"))

# all_dmat_plt
ggsave("all_dmat.png", all_dmat_plt,
       width = 16, height = 16, dpi = 400, scale = .5)

all_edges <- all_dmat %>% 
  select(prs_1, prs_2) %>% 
  mutate(prs_1 = factor(prs_1, unique(prs_1)),
         prs_2 = factor(prs_2, levels(prs_1))) %>% 
  transmute(from = as.integer(prs_1),
            to = as.integer(prs_2))

