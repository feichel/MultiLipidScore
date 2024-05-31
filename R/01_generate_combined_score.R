
# Functions and packages
library(tidyverse)
source(here::here("R/00_functions.R"))


# Estimating effects in DIVAS ---------------------------------------------
divas_raw <- read_rds(here::here("data/Within_class_FA_PE.rds")) |>
  janitor::clean_names() |>
  rename(sex = gender) |>
  mutate(diet = if_else(diet == "A", "SFA", "UFA"))

# Read DIVAS data and log-transform lipid concentrations
divas <- divas_raw |>
  mutate(across(contains("_v"), ~log_transform(., scale = FALSE)))

# Read DIVAS data, log-transform, and scale lipid concentrations
divas_z <- divas_raw |>
  mutate(across(contains("_v"), ~log_transform(., scale = TRUE)))


# Calculate diet effects on log-transformed lipid concentrations and wrangle results data

weights <- divas |>
  select(contains("_v1"),
         -contains("total"))  |>
  names() |>
  tibble(lipid = _) |>
  mutate(lipid = str_remove(lipid, "_v1"),
         model = map(lipid, ~lm.wrapper(divas, .x)),
         coeffs = map(model, broom::tidy)) |>
  unnest(coeffs) |>
  filter(term == "dietUFA") |>
  mutate(p.value_adj = p.adjust(p.value, method = "fdr")) |>
  filter(p.value_adj < 0.05) |>
  divas_lip_names()


# Calculate diet effects on scaled log-transformed lipid concentrations and wrangle results data
weights_z <- divas |>
  select(contains("_v1"),
         -contains("total")) |>
  names() |>
  tibble(lipid = _) |>
  mutate(lipid = str_remove(lipid, "_v1"),
         model = map(lipid, ~lm.wrapper(divas_z, .x)),
         coeffs = map(model, broom::tidy)) |>
  unnest(coeffs) |>
  filter(term == "dietUFA") |>
  mutate(p.value_adj = p.adjust(p.value, method = "fdr")) |>
  filter(p.value_adj < 0.05) |>
  divas_lip_names()


# save weights for scaled and unscaled effects
write_rds(weights, here::here("data/weigths_comb_score.rds"))
write_rds(weights_z, here::here("data/weigths_z_comb_score.rds"))


# Generating score in EPIC ------------------------------------------------

# Read subcohort info
ia_subcohort_id <- feather::read_feather(here::here("data/data_for_analysis.feather")) |>
  select(omics_id, ia_subcohort)


# Read and wrangle EPIC-Potsdam data and scale. Only select lipids that were
# significantly changed by intervention in DIVAS
df_lips <- feather::read_feather("/home/fabian/lipidomics.meta.info/data/ready_dfs/lip_df_imputed.feather") |>
  rename_all(tolower) |>
  rename_all(~str_remove_all(., "^a")) |>
  rename_all(~str_replace_all(., "tag", "tag_")) |>
  rename_all(~str_remove_all(., "_$")) |>
  rename_all(~str_replace_all(., "__", "_")) |>
  select(omics_id, all_of(weights$lipid)) |>
  left_join(ia_subcohort_id) |>
  mutate(across(-c(omics_id, ia_subcohort), ~log_transform(., scale = FALSE)))


df_lips_z <- feather::read_feather("/home/fabian/lipidomics.meta.info/data/ready_dfs/lip_df_imputed.feather") |>
  rename_all(tolower) |>
  rename_all(~str_remove_all(., "^a")) |>
  rename_all(~str_replace_all(., "tag", "tag_")) |>
  rename_all(~str_remove_all(., "_$")) |>
  rename_all(~str_replace_all(., "__", "_")) |>
  select(omics_id, all_of(weights$lipid)) |>
  left_join(ia_subcohort_id) |>
  mutate(across(-c(omics_id, ia_subcohort), ~log_transform(., scale = TRUE)))


# save data for cluster analysis
df_lips |>
  write_rds(here::here("data/data_for_clusters.rds"))

# Get weights from trial
lipids_to_score <- df_lips |>
  select(all_of(weights$lipid))

lipids_to_score_z <- df_lips_z |>
  select(all_of(weights$lipid))

# Calculate score by multiplying weights with respective columns
combined_score <- tibble(omics_id = df_lips$omics_id,
                         comb_score = as.numeric(as.matrix(lipids_to_score) %*% weights$estimate))

combined_score_z <- tibble(omics_id = df_lips$omics_id,
                         comb_score = as.numeric(as.matrix(lipids_to_score_z) %*% weights_z$estimate))



write_rds(combined_score, here::here("data/combined_score.rds"))
write_rds(combined_score_z, here::here("data/combined_score_z.rds"))


# Calculating score in DIVAS ----------------------------------------------
# For estimation of effect on score
weights_renamed <- weights |>
  mutate(lipid = str_replace(lipid, "fa", "fa_") |>
           str_replace("peonly_", "pe_") |>
           str_replace("pe_p", "pep") |>
           str_replace("pe_o", "peo"))

# for baseline
lips_v1 <- divas |>
  select(all_of(paste0(weights_renamed$lipid, "_v1")))

# for post intervention
lips_v2 <- divas |>
  select(all_of(paste0(weights_renamed$lipid, "_v2")))


combined_score_divas <- divas |>
  select(volunteer, diet, age, sex, bmi) |>
  bind_cols(tibble(comb_score_v1 = -as.numeric(as.matrix(lips_v1) %*% weights_renamed$estimate),
                   comb_score_v2 = -as.numeric(as.matrix(lips_v2) %*% weights_renamed$estimate)))

# Save
write_rds(combined_score_divas, here::here("data/combined_score_divas.rds"))



weights_renamed_z <- weights_z |>
  mutate(lipid = str_replace(lipid, "fa", "fa_") |>
           str_replace("peonly_", "pe_") |>
           str_replace("pe_p", "pep") |>
           str_replace("pe_o", "peo"))

# for baseline
lips_v1_z <- divas_z |>
  select(all_of(paste0(weights_renamed$lipid, "_v1")))

# for post intervention
lips_v2_z <- divas_z |>
  select(all_of(paste0(weights_renamed$lipid, "_v2")))


combined_score_divas <- divas_z |>
  select(volunteer, diet, age, sex, bmi) |>
  bind_cols(tibble(comb_score_v1 = -as.numeric(as.matrix(lips_v1_z) %*% weights_renamed_z$estimate),
                   comb_score_v2 = -as.numeric(as.matrix(lips_v2_z) %*% weights_renamed_z$estimate)))

# Save
write_rds(combined_score_divas, here::here("data/combined_score_divas_z.rds"))
