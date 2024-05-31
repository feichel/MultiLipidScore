
# Functions and packages
library(tidyverse)
source(here::here("R/00_functions.R"))



# Load data
diet_scores <- read_rds(here::here("data/separate_scores.rds"))
comb_scores <- read_rds(here::here("data/combined_score.rds"))
comb_scores_z <- read_rds(here::here("data/combined_score_z.rds"))
pheno_df <- feather::read_feather(here::here("data/pheno_df.feather"))



# Combine EPIC data
data_for_analysis <- pheno_df %>%
  # left_join(diet_scores) %>%
  left_join(comb_scores) %>%
  mutate(non_hdl = corr_chol - corr_hdl,
         homa_ir = (corr_insulin * corr_glucose) / 405)


# Combine EPIC data
data_for_analysis_z <- pheno_df %>%
  # left_join(diet_scores) %>%
  left_join(comb_scores_z) %>%
  mutate(non_hdl = corr_chol - corr_hdl,
         homa_ir = (corr_insulin * corr_glucose) / 405)

# save
write_rds(x = data_for_analysis, file = here::here("data/data_for_analysis.rds"))
write_rds(x = data_for_analysis_z, file = here::here("data/data_for_analysis_z.rds"))



