library(tidyverse)


ia_subcohort_id <- read_rds(here::here("data/data_for_analysis.rds")) %>%
  select(omics_id, ia_subcohort)


# read in subcohort lipid concentrations
df_lips <- feather::read_feather("/home/fabian/lipidomics.meta.info/data/ready_dfs/lip_df_imputed.feather") %>%
  rename_all(tolower) %>%
  rename_all(~str_remove_all(., "^a")) %>%
  rename_all(~str_replace_all(., "tag", "tag_")) %>%
  rename_all(~str_remove_all(., "_$")) %>%
  rename_all(~str_replace_all(., "__", "_")) %>%
  left_join(ia_subcohort_id) %>%
  mutate(across(-c(omics_id, ia_subcohort), log)) %>%
  filter(ia_subcohort == 1)

# get list of matched lipids (epic<->predimed) and wrangle names
matching_table <- read_rds(here::here("doc/matching_for_predimed.rds")) %>%
  mutate(across(c(lipid_display, within_class_fa_sum), ~str_to_lower(.) %>%
                  str_replace_all("dg", "dag") %>%
                  str_replace_all("tg", "tag") %>%
                  str_replace_all("\\(", "_") %>%
                  str_remove_all("\\)") %>%
                  str_replace_all(":", "_") %>%
                  str_replace_all("pe_", "peonly_") %>%
                  str_replace_all("pep_", "pe_p_") %>%
                  str_replace_all("peo_", "pe_o_")))

# generate weights (with linear regression within-class fa sum ~ species matches)
results_matching <- matching_table %>%
  filter(!is.na(metabolite)) %>%
  group_by(within_class_fa_sum) %>%
  summarise(matches_for_search = list(lipid_display),
            matches = str_c(lipid_display, collapse = " + ")) %>%
  filter(within_class_fa_sum != matches) %>%
  mutate(formula = str_c(within_class_fa_sum, "~", matches, sep = " "),
         coeffs = map(formula, ~lm(., data = df_lips) %>%
                        broom::tidy() %>%
                        filter(term != "(Intercept)")))

# Read lipid concentrations for all participants
df_lips_for_score <- feather::read_feather("/home/fabian/lipidomics.meta.info/data/ready_dfs/lip_df_imputed.feather") %>%
  rename_all(tolower) %>%
  rename_all(~str_remove_all(., "^a")) %>%
  rename_all(~str_replace_all(., "tag", "tag_")) %>%
  rename_all(~str_remove_all(., "_$")) %>%
  rename_all(~str_replace_all(., "__", "_")) %>%
  left_join(ia_subcohort_id) %>%
  mutate(across(-c(omics_id, ia_subcohort), log))


# function to calculate predicted within-class fa sum (using weights from species)
make_proxy_lip <- function(.variables) {

  select_lips <- results_matching %>%
    select(within_class_fa_sum, coeffs) %>%
    unnest(coeffs) %>%
    filter(within_class_fa_sum == paste(.variables))

  lipids_to_score <- df_lips_for_score %>%
    select(select_lips$term)

  return(as.numeric(as.matrix(lipids_to_score) %*% select_lips$estimate))
}

# Save weights for application in predimed
results_matching %>%
  select(within_class_fa_sum, coeffs) %>%
  unnest(coeffs) %>%
  write_tsv(here::here("data/weights_for_predimed.tsv"))


# calculate predicted within-class fa sum
prox_lips <- map_dfc(results_matching$within_class_fa_sum %>%
                       set_names(.), make_proxy_lip)


# function to assess predicted vs original within-class fa sum correlation
corr_proxy_original <- function(.variables) {
  x <- df_lips_for_score %>%
    select({{.variables}})

  y <- prox_lips %>%
    select({{.variables}})

  as.numeric(cor(x,y))
}

# calculate correlation
corr_res <- map_dfc(
  names(prox_lips) %>%
    set_names(.),
  corr_proxy_original) %>%
  pivot_longer(everything())

# make new data set by exchanging within-class fa sums with predicted within-class fa sums
new_df_lips <- df_lips_for_score %>%
  select(-c(names(prox_lips))) %>%
  bind_cols(prox_lips)

#save to make available for score calculation
write_rds(new_df_lips, here::here("data/lips_for_predimed.rds"))


