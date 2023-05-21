######### OBS - THIS FILE IS ONLY MEANT TO BE RUN USING SOURCE FROM THE RMD DOCUMENT ##############

######## Read and filter segments by type
########
arch_raw <- read_delim("ArchaicSegments.txt", delim = "\t", show_col_types = F) %>% 
  filter(chrom %in% chromosome) %>% 
  filter(length > 0) %>% 
  mutate(country = str_extract(country, "[\\w ]+")) %>% 
  nest(all = !c(region, country)) %>% 
  group_by(country) %>% 
  mutate(
    country = paste0(country, rep("_", n()), region),
    country = if (n() == 1) str_remove(country, "_.+") else country
  ) %>% 
  ungroup %>%
  unnest(all) %>%
  # Prettier region, country and population names
  mutate(
    region = str_replace_all(region, "(?<!^)(?=[A-Z])", "-"),
    across(c(country, pop), ~.x %>% 
             str_replace("_", " (") %>% 
             str_replace("(?<= \\([\\w]{1,25})(?=(\\s|$))", ")")
    )
  ) %>% 
  mutate(chrom = factor(chrom, unique(chrom)[order(as.numeric(unique(chrom)))])) %>% 
  mutate(
    primary_ancestry = apply(across(contains("Shared_with")), 1, which.max),
    primary_ancestry = ifelse(
      rowSums(across(contains("Shared_with"))) < ceiling(snps/2), 
      0, primary_ancestry
    )
  )

arch <- arch_raw %>% 
  filter(primary_ancestry %in% segment_type)

if (epas1) {
  arch_person_sum <- arch %>% 
    distinct(name, pop, country, region, .keep_all = T) %>% 
    mutate(across(!c(name, pop, country, region), ~NA))
  
  arch <- arch %>% 
    filter(chrom == "2" & between(start, 46293667-10^6, 46386697+10^6) & between(end, 46293667-10^6, 46386697+10^6)) %>% 
    nest(dat = !name) %>%
    mutate(name = factor(name, levels = arch_person_sum$name)) %>% 
    complete(name) %>% 
    mutate(dat = map2(dat, name, function(d, n) {
      if (!is.null(d)) return(d)
      
      arch_person_sum[arch_person_sum$name == n, ] %>% 
        select(!name)
    })) %>% 
    unnest(dat) %>% 
    mutate(chrom = "2") %>% 
    ungroup
}

######## Reformat the data structure of the Archaic segments into pseudo-class
######## Approximation is also applied here
all_chrom_range <- arch %>%
  group_by(chrom) %>% 
  summarize(
    start = min(start, na.rm = T),
    end = max(end, na.rm = T),
    .groups = "drop"
  ) %>% 
  arrange(chrom) %>% 
  mutate(
    length = end - start,
    offset = lag(cumsum(length), default = 0) #- first(length)
  )

if (nrow(all_chrom_range) == 1) {
  all_chrom_summary <- as.list(all_chrom_range[-1]) %>% 
    lapply(function(x) set_names(x, all_chrom_range$chrom))
} else {
  all_chrom_summary <- all_chrom_range %>% 
    summarize(across(!chrom, ~list(set_names(.x, chrom)))) %>% 
    as.list %>% 
    lapply(first)
}

all_chrom_range %>% 
  mutate(across(c(start, end), ~.x - start + 1 + offset)) %>% 
  mutate(chrom = as.integer(chrom)) %>%
  mutate(nchrom = lead(chrom),
         nstart = lead(start)) %>% 
  mutate(across(where(is.numeric), ~.x / downscale)) %>% 
  ggplot() +
  geom_segment(aes(chrom, start, xend = chrom, yend = end)) +
  geom_segment(aes(chrom, end, xend = nchrom, yend = nstart)) +
  scale_x_continuous(labels = levels(all_chrom_range$chrom), breaks = 1:length(levels(all_chrom_range$chrom)))


arch_approx <- arch %>% 
  group_by(chrom) %>%
  mutate(
    across(c(start, end), ~
             .x 
           - all_chrom_summary$start[as.character(first(chrom))] 
           + all_chrom_summary$offset[as.character(first(chrom))] 
           + 1)
  ) %>%
  ungroup %>% 
  mutate(
    across(c(start, end), ~.x / downscale),
    start = ceiling(start),
    end = floor(end)
  ) %>%
  mutate(length = end - start) %>% 
  filter(length > 0 | (is.na(start) & is.na(end)))

x_chrom_end_range <- arch_approx %>%
  pull(end) %>% 
  range(na.rm = T)

x_chrom_start_range <- arch_approx %>% 
  pull(start) %>% 
  range(na.rm = T)

x_chrom_min <- x_chrom_start_range[1]
x_chrom_max <- x_chrom_end_range[2] - x_chrom_min + 1

arch_approx <- arch_approx %>%
  mutate(across(c(start, end), ~.x - x_chrom_min + 1))

arch_persons <- arch_approx %>% 
  group_by(name, pop, country, region) %>% 
  arrange(start, end) %>% 
  summarize(
    prs = list(person(start, end, x_chrom_max)),
    .groups = "drop"
  ) %>% 
  mutate(prs = prs %>% 
           set_names(name))

arch_person_length <- arch_persons$prs %>% 
  sapply(
    function(x)
      sum(x$end - x$start)
  )

if (epas1) {
  epas_set <- arch_persons$prs %>% 
    sapply(function(x) sapply(x, function(y) any(is.na(y)))) %>% 
    colSums() %>% 
    equals(0) %>% 
    {names(.)[.]}
} else {
  epas_set <- c()
}

######### Compute pairwise dIoU matrix or load a precomputed one
#########
if (precomputed) {
  all_dmat <- readRDS(paste0("all_dmat", file_suffix, ".rds"))
} else {
  time_used <- system.time(
    all_dmat <- arch_persons %>%
      # slice_sample(n = 25) %>% 
      mutate(prs = trim_non_overlapping_regions(prs)) %>%
      person_pairwise_dist_df("prs", "name", implicit = F)
  )
  
  saveRDS(all_dmat, paste0("all_dmat", file_suffix, ".rds"))
}