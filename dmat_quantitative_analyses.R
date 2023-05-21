############## OBS - ONLY MEANT TO BE RUN AFTER RUNNING THE HYPERPARAMETER CELL IN THE RMD  ########

geo_anova <- all_dmat %>% 
  # filter(region_1 != "Melanesia" & region_2 != "Melanesia") %>% 
  arrange(region_1, region_2, country_1, country_2, pop_1, pop_2, prs_1, prs_2) %>% 
  mutate(
    prs_1 = factor(prs_1),
         prs_2 = factor(prs_2, levels(prs_1))
    ) %>% 
  mutate(
    pop_1     = ifelse(pop_1 == pop_2, pop_1, "Between") %>% 
      factor(c("Between", sort(unique(pop_1)))),
    
    country_1 = ifelse(country_1 == country_2, country_1, "Between") %>% 
      factor(c("Between", sort(unique(country_1)))),
    
    region_1  = ifelse(region_1 == region_2, region_1, "Between") %>% 
      factor(c("Between", sort(unique(region_1)))),
    
    pop_2     = ifelse(pop_2 == pop_1, pop_2, "Between") %>% 
      factor(c("Between", sort(unique(pop_2)))),
    
    country_2 = ifelse(country_2 == country_1, country_2, "Between") %>% 
      factor(c("Between", sort(unique(country_2)))),
    
    region_2  = ifelse(region_2 == region_1, region_2, "Between") %>% 
      factor(c("Between", sort(unique(region_2))))
  ) %>% 
  filter(as.integer(prs_1) < as.integer(prs_2)) %>% 
  select(iou, contains("_1")) %>%
  # group_by(pop_1) %>% 
  # filter(length(unique(prs_1)) > 1 | first(pop_1) == "Between") %>% 
  # group_by(country_1) %>% 
  # filter(length(unique(pop_1)) > 1 | first(country_1) == "Between") %>% 
  # group_by(region_1) %>% 
  # filter(length(unique(country_1)) > 1 | first(region_1) == "Between") %>% 
  # ungroup %>% 
  rename_with(function(x) str_remove(x, "_\\d$")) %>%
  mutate(between = rowSums(across(c(region, country, pop), ~ .x == "Between")) %>% 
           factor) %>% 
  # mutate(iou = log(iou/(1 - iou))) %>%
  lm(iou ~ I(region == "Between") + region
         + I(country == "Between") + country
         + I(pop == "Between") + pop, 
     data = .)

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

all_dmat %>% 
  filter(!epas1 | (prs_1 %in% epas_set & prs_2 %in% epas_set)) %>% 
  mutate(
    len_1 = arch_person_length[prs_1] * downscale,
    len_2 = arch_person_length[prs_2] * downscale
  ) %>% 
  mutate(region_1 = if (!epas1) region_1 else ifelse(country_1 == "Tibet", "Tibet", region_1),
         region_2 = if (!epas1) region_2 else ifelse(country_2 == "Tibet", "Tibet", region_2)) %>%
  filter(region_1 == region_2) %>%
  arrange(region_1, region_2, country_1, country_2, pop_1, pop_2, prs_1, prs_2) %>%
  mutate(prs_1 = factor(prs_1, unique(prs_1)),
         prs_2 = factor(prs_2, levels(prs_1))) %>%
  filter(as.numeric(prs_1) > as.numeric(prs_2)) %>%
  mutate(
    len = len_1,
    # len_1 = cut(len_1, breaks = scales::breaks_pretty(n = 30)(range(len)), dig.lab = 25),
    # len_2 = cut(len_2, breaks = scales::breaks_pretty(n = 30)(range(len)), dig.lab = 25)
    len_1 = cut_width(len_1, {if (!epas1) 10000 else 50000} * downscale, boundary = 0, dig.lab = 25),
    len_2 = cut_width(len_2, {if (!epas1) 10000 else 50000} * downscale, boundary = 0, dig.lab = 25)
  ) %>%
  filter(len_1 == len_2) %>%
  group_by(len_1, region_1) %>%
  summarize(
    iou = median(iou),
    n = n(),
    .groups = "drop"
  ) %>% 
  mutate(len_1 = droplevels(len_1)) %>%  
  complete(len_1, region_1, fill = list(n = 0)) %>% 
  mutate(
    iou = ifelse(is.na(iou), "-", scales::label_percent(.1)(iou)),
    iou_n = paste0(iou, " (\\textit{", n, "})"),
    region_1 = region_1 %>% 
      factor(c("West-Eurasia", 
               "South-Asia", 
               "East-Asia", 
               "Tibet", 
               "Central-Asia-Siberia", 
               "Melanesia")),
    region_1 = droplevels(region_1)
  ) %>% 
  arrange(region_1) %>% 
  select(!c(iou,n)) %>% 
  mutate(
    iou_n = str_replace(iou_n, "%", "\\\\%"),
    len_1 = if (!epas1) str_replace_all(len_1, "0{6}(?!0)", "M") else str_replace_all(len_1, "0{3}(?!0)", "K"),
    len_1 = len_1 %>% 
      str_replace("(?<=(\\(|\\[)\\d{2}[MK],)(?=\\d{3}[MK])", "\\\\hspace{0.25cm}") %>% 
      str_replace("(?<=(\\(|\\[)\\d{2}[MK],)(?=\\d{2}[MK])", "\\\\hspace{0.5cm}") %>% 
      str_replace("(?<=(\\(|\\[)0,)", "\\\\hspace{0.75cm}") %>% 
      str_replace("(?<=(\\(|\\[)\\d{3}[MK],)", " "),
  ) %>% 
  mutate(
    region_1 = paste0("\\textbf{", region_1, "}")
  ) %>% 
  pivot_wider(
    id_cols = len_1,
    names_from = region_1,
    values_from = iou_n
  ) %>% 
  rename("\\makecell{\\textbf{Total Archaic}\\\\\\textbf{Segment Length}}" = 1) %>% 
  kable("latex",
        booktabs = T,
        escape = F,
        align = "c") %>%
  add_header_above(c("", "$Q_{\\\\%50}(d_\\\\mathrm{IoU})$  (\textit{n})" = if (!epas1) 5 else 6), escape = F) %>% 
  as.character %>% 
  str_split("\\n") %>% 
  unlist %>% 
  extract(which(nchar(.) != 0)) %>% 
  paste0(collapse = "\n\t") %>%
  {paste0("\t", .)} %>% 
  write_lines(paste0("segment_length_cor", file_suffix,".txt"))
  


# all_coef <- all_dmat %>%
#   filter(!epas1 | (prs_1 %in% epas_set & prs_2 %in% epas_set)) %>% 
#   arrange(region_1, region_2, country_1, country_2, pop_1, pop_2, prs_1, prs_2) %>%
#   mutate(prs_1 = factor(prs_1),
#          prs_2 = factor(prs_2, levels(prs_1)),) %>%
#   mutate(
#     len_1 = arch_person_length[prs_1] * 1000,
#     len_2 = arch_person_length[prs_2] * 1000
#   ) %>% 
#   nest(dat = !contains("_1")) %>% 
#   mutate(
#     mod1 = map(dat, function(x) lmer(iou ~ len_2 + (1 | region_2) + (1 | country_2) + (1 | pop_2), data = x))
#   ) %>% 
#   mutate(
#     coef = map(mod1, function(x) as_tibble(set_colnames(t(coefficients(summary(x))[1:2,1]), c("Intercept", "Slope"))))
#   ) %>% 
#   unnest(coef) %>% 
#   select(!c(dat, mod1)) 
# 
# tibble(
#   coef = list(all_coef %>% 
#                 lmer(Slope ~ len_1 + (1 | region_1 / country_1 / pop_1), data = .) %>% 
#                 summary %>% 
#                 coefficients() %>% 
#                 as.data.frame %>% 
#                 rownames_to_column("second_term") %>% 
#                 as_tibble,
#               all_coef %>% 
#                 lmer(Intercept ~ len_1 + (1 | region_1 / country_1 / pop_1), data = .) %>% 
#                 summary %>% 
#                 coefficients() %>% 
#                 as.data.frame %>% 
#                 rownames_to_column("second_term") %>% 
#                 as_tibble
#   ),
#   term = c("Intercept", "Slope")
# ) %>% 
#   unnest(coef) %>% 
#   relocate(term, .before = second_term) %>% 
#   mutate(second_term = c("(Intercept)" = "Intercept",
#                          "len_1" = "Slope")[second_term]) %>% 
#   mutate(across(where(is.numeric), label_prettify_scientific(F,3)))
# 
# 
# all_coef %>% 
#   pivot_longer(!contains("_1"), names_to = "parm", values_to = "est") %>% 
#   ggplot(aes(len_1, est, color = region_1)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   scale_y_continuous(trans = abslog(),
#                      labels = prettify_scientific) +
#   scale_x_continuous(labels = prettify_scientific) +
#   facet_wrap(~parm, scales = "free_y")
