manual_2d_density <- arch %>% 
  mutate(
    across(contains("Shared"), ~.x/snps)
  ) %>% 
  select(c(region, contains("Shared"))) %>% 
  filter(rowSums(across(!region) == 0) == 0) %>%
  group_by(region) %>% 
  reframe(
    density = MASS::kde2d(Shared_with_Altai, Shared_with_Denisova, n = 100)
  )


manual_2d_density %>% 
  mutate(
    col_name = names(density)
  ) %>% 
  pivot_wider(
    id_cols = region,
    names_from = col_name,
    values_from = density
  ) %>% 
  mutate(
    density = pmap(list(x,y,z), function(xc, yc, d) {
      tibble(
        density = as.vector(d),
        x = rep(xc, length(xc)),
        y = rep(yc, each = length(yc))
      )
    })
  ) %>% 
  select(region, density) %>% 
  unnest(density) %>% 
  mutate(
    density = ifelse(density < 0.01, NA, density)
  ) %>% 
  ggplot(aes(x,y,fill=density)) +
  geom_raster() +
  scale_fill_viridis_c(trans = "log10") +
  coord_equal(expand = F) +
  facet_wrap(~region)
