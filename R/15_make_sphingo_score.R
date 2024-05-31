# to do get MLS-ceramides that overlap with lipogain-2
# generate score
# assess divas effect on score
# --> scale in epic by divas observed effect size
# assess diet correlations in epic with cer-score
# assess risk assoc with in epic with cer-score

library(tidyverse)

divas <- read_rds(here::here("data/Within_class_FA_PE.rds")) |>
  janitor::clean_names() |>
  rename(sex = gender) |>
  mutate(diet = if_else(diet == "A", "SFA", "UFA"))

lipogain_2 <- c("sm_fa14_0", "cer_fa18_0", "dcer_fa18_0",
                "dcer_fa20_0", "dcer_fa22_0", "dcer_fa24_1", "hcer_fa18_0")



weights <- read_rds(here::here("data/weigths_comb_score.rds")) |>
  filter(lipid %in% lipogain_2) |>
  mutate(lipid = str_replace(lipid, "fa", "fa_"))


weights_z <- read_rds(here::here("data/weigths_z_comb_score.rds")) |>
  filter(lipid %in% lipogain_2) |>
  mutate(lipid = str_replace(lipid, "fa", "fa_"))


# for baseline
lips_v1 <- divas |>
  select(all_of(paste0(weights$lipid, "_v1")))

# for post intervention
lips_v2 <- divas |>
  select(all_of(paste0(weights$lipid, "_v2")))


sphingolipid_score <- divas |>
  select(volunteer, diet, age, sex, bmi) |>
  bind_cols(tibble(sl_score_v1 = -as.numeric(as.matrix(lips_v1) %*% weights$estimate),
                   sl_score_v2 = -as.numeric(as.matrix(lips_v2) %*% weights$estimate)))

# Save


# Calculate effects on score and clinical markers
divas_effects <- tibble(var = "sl_score") %>%
  mutate(model = map(var, ~lm.wrapper(sphingolipid_score, .)),
         coeffs = map(model, ~broom::tidy(., conf.int = TRUE, conf.level = 0.95))) %>%
  unnest(coeffs) %>%
  filter(str_detect(term, "diet")) %>%
  mutate(p.value_adj = p.adjust(p.value, method = "fdr"))

# sl_score weight  = -0.355



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
  select(omics_id, all_of(str_replace(weights$lipid, "fa_", "fa"))) |>
  left_join(ia_subcohort_id) |>
  mutate(across(-c(omics_id, ia_subcohort), ~log_transform(., scale = FALSE)))


# Get weights from trial
lipids_to_score <- df_lips |>
  select(all_of(str_replace(weights$lipid, "fa_", "fa")))



# Calculate score by multiplying weights with respective columns
sphingolipid_score_epic <- tibble(omics_id = df_lips$omics_id,
                                  sphingolipid_score = as.numeric(as.matrix(lipids_to_score) %*% weights$estimate))


# Load data
pheno_df <-  read_rds(here::here("data/data_for_analysis.rds")) %>%
  mutate(across(c(epic_cvd, case_diab_caco, fasting, educc3, smk_4cat, alccat,
                  antihyp, lipidlower, ass), as_factor)) %>%
  left_join(sphingolipid_score_epic) |>
  mutate(across(c(
    "corr_trigly",
    "corr_chol",
    "corr_hdl",
    "corr_hscrp",
    "non_hdl"), ~log(.) %>%
      #scale() %>%
      as.numeric()),
    sphingolipid_score = sphingolipid_score/0.355) # this is the effect seen in DIVAS



write_rds(pheno_df, here::here("data/sphingolipid_score_epic.rds"))

