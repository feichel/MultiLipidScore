library(tidyverse)



df_figure_weights <- read_rds(here::here("data/lipids_and_weights.rds"))



fa_infos <-  xlsx::read.xlsx(here::here("doc/fa_names_list.xlsx"), sheetIndex = 1) %>%
  tibble() %>%
  janitor::clean_names() %>%
  mutate(shorthand = str_remove(shorthand, "C") %>%
           str_squish(),
         trans_or_cis = if_else(is.na(trans_or_cis), "saturated", trans_or_cis),
         common_name = if_else(trans_or_cis == "saturated", common_name,  str_c(common_name, " (", omega_n, ")"))) %>%
  filter(trans_or_cis != "trans") %>%
  select(shorthand, common_name) %>%
  group_by(shorthand) %>%
  summarize(common_name = paste(common_name, collapse = ", "))


df_figure_weights %>%
  select(lipid = lipid_display, weight) %>%
  separate(lipid, remove = FALSE, into = c("class", "l", "db"),
           extra = "drop") %>%
  mutate(shorthand = str_c(l, ":", db) %>%
           str_remove("FA"),
         unit = "(µmol/L)") %>%
  left_join(fa_infos) %>%
  select(lipid, class, fa_short_hand = shorthand, common_name, weight, unit) %>%
  xlsx::write.xlsx(., "/home/fabian/alle_shortcut/!MEP/Projekte/EPIC-Potsdam/Diabetes/Lipidomics/MultiLipidScore/Revision_1/supplemental_table_7_8.xlsx")


# full formula for copy and paste
df_figure_weights |>
  select(lipid_display, estimate) |>
  mutate(estimate = str_c("(", round(estimate, 2), ")", "*", lipid_display, " (µmol/L)")) |>
  pivot_longer(-lipid_display) |>
  summarize(formula = str_c(value, collapse = " + "), .by = name) |>
  pull(formula)
