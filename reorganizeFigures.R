library(tidyverse)


tibble(
  file = list.files()
) %>% 
  mutate(
    file_type = str_extract(file, "(?<=\\.)[a-zA-Z]{2,4}$") %>% 
      tolower(),
    file_type = ifelse(nchar(file) == 0, "directory", file_type)
  ) %>% 
  filter(file_type == "pdf") %>% 
  mutate(
    figure_type = str_extract(tolower(file), "(?<=_)[a-zA-Z0-9]+(?=\\.(pdf|png)$)"),
    figure_type = ifelse(str_detect(file, "^approx"), "approx", figure_type),
    figure_type = ifelse(figure_type %in% c("altai", "denisova"   , "vindija", 
                                            "known", "neanderthal", "unknown", 
                                            "altaidenisova", "epas1", "approx"),
                         figure_type, "all")
  ) %>% 
  select(!file_type) %>% 
  nest(files = !figure_type) %>% 
  mutate(files = map(files, unlist)) %>% 
  mutate(
    temp = map2(figure_type, files, function(d, f) {
      if (!file.exists("figures")) dir.create("figures")
      if (!file.exists(paste0("figures/", d))) dir.create(paste0("figures/", d))
      
      file.rename(
        from = paste0(f),
        to = paste0("figures/", d, "/", f)
      )
    })
  )
