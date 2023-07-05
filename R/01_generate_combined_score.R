
# Functions and packages
library(tidyverse)
source(here::here("R/00_functions.R"))


# Estimating effects in DIVAS ---------------------------------------------

# Read DIVAS data and logtransform lipid concentrations
divas <- read_rds(here::here("data/Within_class_FA_PE.rds")) %>%
  janitor::clean_names() %>%
  mutate(diet = if_else(diet == "A", "SFA", "UFA"),
         across(contains("_v"), ~log(.) %>%
                  #scale() %>%
                  as.numeric())) %>%
  rename(sex = gender)


# Read DIVAS data, log-transform, and scale lipid concentrations
divas_z <- read_rds(here::here("data/Within_class_FA_PE.rds")) %>%
  janitor::clean_names() %>%
  mutate(diet = if_else(diet == "A", "SFA", "UFA"),
         across(contains("_v"), ~log(.) %>%
                  scale() %>%
                  as.numeric())) %>%
  rename(sex = gender)

# Calculate diet effects on log-transformed lipid concentrations and wrangle results data
weights <- divas %>%
  select(contains("_v1"),
         -contains("total")) %>%
  names() %>%
  tibble(lipid = .) %>%
  mutate(lipid = str_remove(lipid, "_v1"),
         model = map(lipid, ~lm.wrapper(divas, .x)),
         coeffs = map(model, broom::tidy)) %>%
  unnest(coeffs) %>%
  filter(term == "dietUFA") %>%
  mutate(p.value_adj = p.adjust(p.value, method = "fdr")) %>%
  filter(p.value_adj < 0.05) %>%
  separate(lipid, into = c("base", "trash", "length", "db"), remove = FALSE) %>%
  mutate(lipid_display = if_else(base %in% c("tag", "dag", "pe", "pep", "peo",
                                             "pi", "pc"),
                                 str_c(base, "(FA", length, ":", db, ")"),
                                 str_c(base, "(", length, ":", db, ")")) %>%
           str_to_upper(),
         lipid_display = str_replace(lipid_display, "TAG", "TG") %>%
           str_replace("MAG", "MG") %>%
           str_replace("DAG", "DG") %>%
           str_replace("LCER", "LacCER") %>%
           str_replace("HCER", "HexCER") %>%
           str_replace("DCER", "dhCER") %>%
           str_replace("CER", "Cer"),
         lipid = str_replace(lipid, "fa_", "fa") %>%
           str_replace("pe_fa", "peonly_fa") %>%
           str_replace("peo_", "pe_o_") %>%
           str_replace("pep_", "pe_p_"))


# Calculate diet effects on scaled log-transformed lipid concentrations and wrangle results data
weights_z <- divas %>%
  select(contains("_v1"),
         -contains("total")) %>%
  names() %>%
  tibble(lipid = .) %>%
  mutate(lipid = str_remove(lipid, "_v1"),
         model = map(lipid, ~lm.wrapper(divas_z, .x)),
         coeffs = map(model, broom::tidy)) %>%
  unnest(coeffs) %>%
  filter(term == "dietUFA") %>%
  mutate(p.value_adj = p.adjust(p.value, method = "fdr")) %>%
  filter(p.value_adj < 0.05) %>%
  separate(lipid, into = c("base", "trash", "length", "db"), remove = FALSE) %>%
  mutate(lipid_display = if_else(base %in% c("tag", "dag", "pe", "pep", "peo",
                                             "pi", "pc"),
                                 str_c(base, "(FA", length, ":", db, ")"),
                                 str_c(base, "(", length, ":", db, ")")) %>%
           str_to_upper(),
         lipid_display = str_replace(lipid_display, "TAG", "TG") %>%
           str_replace("MAG", "MG") %>%
           str_replace("DAG", "DG") %>%
           str_replace("LCER", "LacCER") %>%
           str_replace("HCER", "HexCER") %>%
           str_replace("DCER", "dhCER") %>%
           str_replace("CER", "Cer"),
         lipid = str_replace(lipid, "fa_", "fa") %>%
           str_replace("pe_fa", "peonly_fa") %>%
           str_replace("peo_", "pe_o_") %>%
           str_replace("pep_", "pe_p_"))



# save weights for scaled and unscaled effects
write_rds(weights, here::here("data/weigths_comb_score.rds"))
write_rds(weights_z, here::here("data/weigths_z_comb_score.rds"))



# Generating score in EPIC ------------------------------------------------

# Read subcohort info
ia_subcohort_id <- feather::read_feather(here::here("data/data_for_analysis.feather")) %>%
  select(omics_id, ia_subcohort)


# Read and wrangle EPIC-Potsdam data and scale. Only select lipids that were
# significantly changed by intervention in DIVAS
df_lips <- feather::read_feather("/home/fabian/lipidomics.meta.info/data/ready_dfs/lip_df_imputed.feather") %>%
  rename_all(tolower) %>%
  rename_all(~str_remove_all(., "^a")) %>%
  rename_all(~str_replace_all(., "tag", "tag_")) %>%
  rename_all(~str_remove_all(., "_$")) %>%
  rename_all(~str_replace_all(., "__", "_")) %>%
  select(omics_id, all_of(weights$lipid)) %>%
  left_join(ia_subcohort_id) %>%
  mutate(across(-c(omics_id, ia_subcohort), log))


df_lips %>%
  write_rds(here::here("data/data_for_clusters.rds"))

# Get weights from trial (unadjusted for class sums)
lipids_to_score <- df_lips %>%
  select(all_of(weights$lipid))


# Calculate score by multiplying weights with respective columns
combined_score <- tibble(omics_id = df_lips$omics_id,
                         comb_score = as.numeric(as.matrix(lipids_to_score) %*% weights$estimate))


write_rds(combined_score, here::here("data/combined_score.rds"))




# Calculating score in DIVAS ----------------------------------------------

weights_renamed <- weights %>%
  mutate(lipid = str_replace(lipid, "fa", "fa_") %>%
           str_replace("peonly_", "pe_") %>%
           str_replace("pe_p", "pep") %>%
           str_replace("pe_o", "peo"))

# for baseline
lips_v1 <- divas %>%
  select(all_of(paste0(weights_renamed$lipid, "_v1")))

# for post intervention
lips_v2 <- divas %>%
  select(all_of(paste0(weights_renamed$lipid, "_v2")))


combined_score_divas <- divas %>%
  select(volunteer, diet, age, sex, bmi) %>%
  bind_cols(tibble(comb_score_v1 = -as.numeric(as.matrix(lips_v1) %*% weights_renamed$estimate),
                   comb_score_v2 = -as.numeric(as.matrix(lips_v2) %*% weights_renamed$estimate)))

# Save
write_rds(combined_score_divas, here::here("data/combined_score_divas.rds"))


