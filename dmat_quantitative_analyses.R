all_dmat_dist <- all_dmat %>% 
  select(1:3) %>% 
  arrange(prs_1) %>% 
  pivot_wider(
    id_cols = prs_1, names_from = prs_2, values_from = iou
  ) %>% 
  column_to_rownames("prs_1") %>% 
  as.dist


# all_dmat_geocat <- all_dmat %>% 
#   filter(prs_2 == first(prs_2)) %>% 
#   select(prs_1, pop_1, country_1, region_1) %>%
#   arrange(region_1, country_1, pop_1) %>% 
#   mutate(across(!prs_1, ~factor(.x, unique(.x)))) %>% 
#   arrange(prs_1) %>% 
#   select(!prs_1) %>% 
#   rename_with(function(x) str_remove(x, "_\\d$")) %>% 
#   with(interaction(region, country))
# 
# 
# adon_geo <- adonis2(
#   all_dmat_dist ~ all_dmat_geocat,
#   permutations = 10,
#   # strata = all_dmat_geocat
# )
# 
# anosim_geo <- anosim(
#   all_dmat_dist,
#   all_dmat_geocat,
#   permutations = 1000
# )



geo_anova <- all_dmat %>% 
  # filter(region_1 != "Melanesia" & region_2 != "Melanesia") %>% 
  arrange(region_1, region_2, country_1, country_2, pop_1, pop_2, prs_1, prs_2) %>% 
  mutate(prs_1 = factor(prs_1),
         prs_2 = factor(prs_2, levels(prs_1)),) %>% 
  mutate(
    pop_1     = factor(ifelse(pop_1     == pop_2,     pop_1,     "Between"), c("Between", sort(unique(pop_1)))),
    country_1 = factor(ifelse(country_1 == country_2, country_1, "Between"), c("Between", sort(unique(country_1)))),
    region_1  = factor(ifelse(region_1  == region_2,  region_1,  "Between"), c("Between", sort(unique(region_1)))),
    pop_2     = factor(ifelse(pop_2     == pop_1,     pop_2,     "Between"), c("Between", sort(unique(pop_2)))),
    country_2 = factor(ifelse(country_2 == country_1, country_2, "Between"), c("Between", sort(unique(country_2)))),
    region_2  = factor(ifelse(region_2  == region_1,  region_2,  "Between"), c("Between", sort(unique(region_2))))
  ) %>% 
  filter(as.integer(prs_1) < as.integer(prs_2)) %>% 
  select(iou, contains("_1")) %>%
  group_by(pop_1) %>% 
  filter(length(unique(prs_1)) > 1 | first(pop_1) == "Between") %>% 
  group_by(country_1) %>% 
  filter(length(unique(pop_1)) > 1 | first(country_1) == "Between") %>% 
  group_by(region_1) %>% 
  filter(length(unique(country_1)) > 1 | first(region_1) == "Between") %>% 
  ungroup %>% 
  rename_with(function(x) str_remove(x, "_\\d$")) %>%
  mutate(between = rowSums(across(c(region, country, pop), ~ .x == "Between")) %>% 
           factor) %>% 
  # mutate(iou = log(iou/(1 - iou))) %>%
  lm(iou ~ I(region == "Between") + region + I(country == "Between") + country + I(pop == "Between") + pop, data = .)

geo_anova %>% 
  residuals %>% 
  hist(1000, plot = F) %>% 
  unclass %>% 
  lapply(list) %>% 
  as_tibble %>% 
  {.[,sapply(., function(x) length(x[[1]])) == length(.$mids[[1]])]} %>% 
  unnest(everything()) %>%
  filter(between(mids, -.1, .2)) %>%
  ggplot(aes(mids, counts)) +
  geom_col(color = "gray35") +
  coord_cartesian(expand = F)

geo_anova %>% 
  anova(test = "F") %>% 
  broom::tidy() %>%
  mutate(
    pdev = scales::label_percent(.01)(sumsq / sum(sumsq)),
    meansq = meansq / sum(meansq)
  ) %>% 
  mutate(
    between = str_detect(term, "^I\\("),
    term = ifelse(between, str_extract(term, "(?<=^I\\()[a-zA-Z]+"), term),
    between = ifelse(between, "between", "within"),
    pdev = set_names(pdev, term)
  ) %>% 
  select(between, term, pdev) %>% 
  complete(term, between) %>% 
  pivot_wider(
    id_cols = term, names_from = between, values_from = pdev
  ) %>% 
  mutate(term = c("region" = "Region",
                  "country" = "Country",
                  "pop" = "Population",
                  "Residuals" = "Residuals")[term],
         term = factor(term, c("Region", "Country", "Population", "Residuals"))) %>% 
  arrange(term) %>% 
  kable(
    "latex",
    booktabs = T,
    col.names = c("", "Between", "Within"),
    escape = F
  ) %>% 
  row_spec(0, bold = T) %>% 
  add_header_above(c("", "Proportion of $d_\\mathrm{IoU}$ explained by" = 2), escape = F, bold = T) %>% 
  write_lines(paste0("shared_ancestry_explained_geography", file_suffix, ".txt"))
 # Average IoU for individuals both with total segment lengths within a bin

arch_person_length <- arch_persons$prs %>% 
  sapply(
    function(x)
      sum(x$end - x$start)
  )

all_dmat %>%
  filter(!epas1 | (prs_1 %in% epas_set & prs_2 %in% epas_set)) %>% 
  arrange(region_1, region_2, country_1, country_2, pop_1, pop_2, prs_1, prs_2) %>%
  mutate(prs_1 = factor(prs_1),
         prs_2 = factor(prs_2, levels(prs_1))) %>%
  mutate(
    len_1 = arch_person_length[prs_1] * downscale,
    len_2 = arch_person_length[prs_2] * downscale
  ) %>%
  mutate(
    len_1 = cut_width(len_1, {if (!epas1) 25000 else 50000} * downscale, boundary = 0, dig.lab = 25),
    len_2 = cut_width(len_2, {if (!epas1) 25000 else 50000} * downscale, boundary = 0, dig.lab = 25)
  ) %>%
  filter(len_1 == len_2) %>%
  filter(as.numeric(prs_1) > as.numeric(prs_2)) %>% 
  group_by(len_1, region_1) %>%
  summarize(
    iou = median(iou),
    n = n(),
    .groups = "drop"
  ) %>% 
  mutate(len_1 = droplevels(len_1)) %>%  
  complete(len_1, region_1, fill = list(n = 0)) %>% 
  mutate(iou = ifelse(is.na(iou), "-", scales::label_percent(.1)(iou)),
         iou_n = paste0(iou, " (\\textit{", n, "})"),
         region_1 = factor(region_1, c("West-Eurasia", "South-Asia", "East-Asia", "Central-Asia-Siberia", "Melanesia"))) %>% 
  arrange(region_1) %>% 
  select(!c(iou,n)) %>% 
  mutate(iou_n = str_replace(iou_n, "%", "\\\\%"),
         len_1 = if (!epas1) str_replace_all(len_1, "0{6}(?!0)", "M") else str_replace_all(len_1, "0{3}(?!0)", "K")) %>% 
  pivot_wider(
    id_cols = len_1,
    names_from = region_1,
    values_from = iou_n
  ) %>% 
  rename("Total Archaic Segment Length" = 1) %>% 
  kable("latex",
        booktabs = T,
        escape = F,
        align = "c") %>%
  row_spec(0, bold = T) %>% 
  add_header_above(c("", "$Q_{\\\\%50}(d_\\\\mathrm{IoU})$" = 5), escape = F) %>% 
  write_lines(paste0("segment_length_cor", file_suffix,".txt"))
  



all_coef <- all_dmat %>%
  arrange(region_1, region_2, country_1, country_2, pop_1, pop_2, prs_1, prs_2) %>%
  mutate(prs_1 = factor(prs_1),
         prs_2 = factor(prs_2, levels(prs_1)),) %>%
  mutate(
    len_1 = arch_person_length[prs_1] * 1000,
    len_2 = arch_person_length[prs_2] * 1000
  ) %>% 
  nest(dat = !contains("_1")) %>% 
  mutate(
    mod1 = map(dat, function(x) lmer(iou ~ len_2 + (1 | region_2) + (1 | country_2) + (1 | pop_2), data = x))
  ) %>% 
  mutate(
    coef = map(mod1, function(x) as_tibble(set_colnames(t(coefficients(summary(x))[1:2,1]), c("Intercept", "Slope"))))
  ) %>% 
  unnest(coef) %>% 
  select(!c(dat, mod1)) 

tibble(
  coef = list(all_coef %>% 
                lmer(Slope ~ len_1 + (1 | region_1 / country_1 / pop_1), data = .) %>% 
                summary %>% 
                coefficients() %>% 
                as.data.frame %>% 
                rownames_to_column("second_term") %>% 
                as_tibble,
              all_coef %>% 
                lmer(Intercept ~ len_1 + (1 | region_1 / country_1 / pop_1), data = .) %>% 
                summary %>% 
                coefficients() %>% 
                as.data.frame %>% 
                rownames_to_column("second_term") %>% 
                as_tibble
  ),
  term = c("Intercept", "Slope")
) %>% 
  unnest(coef) %>% 
  relocate(term, .before = second_term) %>% 
  mutate(second_term = c("(Intercept)" = "Intercept",
                         "len_1" = "Slope")[second_term]) %>% 
  mutate(across(where(is.numeric), label_prettify_scientific(F,3)))


all_coef %>% 
  pivot_longer(!contains("_1"), names_to = "parm", values_to = "est") %>% 
  ggplot(aes(len_1, est, color = region_1)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_y_continuous(trans = abslog(),
                     labels = prettify_scientific) +
  scale_x_continuous(labels = prettify_scientific) +
  facet_wrap(~parm, scales = "free_y")
