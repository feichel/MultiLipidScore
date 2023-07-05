
# Functions and packages
library(tidyverse)
source(here::here("R/00_functions.R"))



# Load data
diet_scores <- read_rds(here::here("data/separate_scores.rds"))
comb_scores <- read_rds(here::here("data/combined_score.rds"))
pheno_df <- feather::read_feather(here::here("data/pheno_df.feather"))

# Combine EPIC data
data_for_analysis <- pheno_df %>%
  left_join(diet_scores) %>%
  left_join(comb_scores) %>%
  mutate(non_hdl = corr_chol - corr_hdl,
         homa_ir = (corr_insulin * corr_glucose) / 405)


# Not used
# data_for_analysis <- data_for_analysis_temp %>%
#   mutate(non_hdl_mufa = (scale(corr_chol - corr_hdl) %>%
#                            as.numeric()) / -0.359,
#          non_hdl_pufa = (scale(corr_chol - corr_hdl) %>%
#                            as.numeric()) / -0.445,
#          mufa_score_divas = mufa_score/29,
#          pufa_score_divas = pufa_score/10,
#          across(c("mufa_score", "pufa_score"), scale_subcohort,
#                 .names = "{.col}_z"))

# save
write_rds(data_for_analysis, here::here("data/data_for_analysis.rds"))


