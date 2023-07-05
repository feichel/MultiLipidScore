library(tidyverse)
source(here::here("R/00_functions.R"))


# Read DIVAS data
divas <- read_rds(here::here("data/Within_class_FA_PE.rds")) %>%
  janitor::clean_names() %>%
  mutate(across(contains("_v"), ~log(.) %>%
                  # scale() %>%
                  as.numeric())) %>%
  rename(sex = gender)


# Calculate effects and wrangle data
weights <- divas %>%
  select(contains("_v1"),
         -contains("total")) %>%
  names() %>%
  tibble(lipid = .) %>%
  mutate(lipid = str_remove(lipid, "_v1"),
         model = map(lipid, ~lm.wrapper(divas, .x)),
         coeffs = map(model, broom::tidy)) %>%
  unnest(coeffs) %>%
  mutate(p.value_adj = p.adjust(p.value, method = "fdr")) %>%
  filter(str_detect(term, "diet"),
         p.value_adj < 0.05) %>%
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
           str_replace("pep_", "pe_p_"),
         term = if_else(term == "dietC", "Mixed UFA", "MUFA"))



# save weights
write_rds(weights, here::here("data/weigths_separate_scores.rds"))


# Make intervention-specific data frames
mufa_weights <- weights %>%
  filter(term == "MUFA") %>%
  na.omit()

pufa_weights <- weights %>%
filter(term == "Mixed UFA") %>%
  na.omit()



# Read subcohort info
ia_subcohort_id <- feather::read_feather(here::here("data/data_for_analysis.feather")) %>%
  select(omics_id, ia_subcohort)


# Read and wrangle EPIC-Potsdam data and scale
df_lips <- feather::read_feather("/home/fabian/lipidomics.meta.info/data/ready_dfs/lip_df_imputed.feather") %>%
  rename_all(tolower) %>%
  rename_all(~str_remove_all(., "^a")) %>%
  rename_all(~str_replace_all(., "tag", "tag_")) %>%
  rename_all(~str_remove_all(., "_$")) %>%
  rename_all(~str_replace_all(., "__", "_")) %>%
  select(omics_id, all_of(weights$lipid)) %>%
  left_join(ia_subcohort_id) #%>%
  # mutate(across(-omics_id, ~scale_subcohort(cur_data(), .)))


# Get weights from trial
mufa_lipids <- df_lips %>%
  select(all_of(mufa_weights$lipid))

pufa_lipids <- df_lips %>%
  select(all_of(pufa_weights$lipid))


# Calculate score by multiplying weights with respective columns
diet_scores <- tibble(omics_id = df_lips$omics_id,
                      mufa_score = as.numeric(as.matrix(mufa_lipids) %*% mufa_weights$estimate),
                      pufa_score = as.numeric(as.matrix(pufa_lipids) %*% pufa_weights$estimate))



# Read EPIC-Potsdam lipidomics (only imputed, not log-transformed or scaled)
df_lips <- feather::read_feather("/home/fabian/lipidomics.meta.info/data/ready_dfs/lip_df_imputed.feather") %>%
  rename_all(tolower) %>%
  rename_all(~str_remove_all(., "^a")) %>%
  rename_all(~str_replace_all(., "tag", "tag_")) %>%
  rename_all(~str_remove_all(., "_$")) %>%
  rename_all(~str_replace_all(., "__", "_"))



# Get info on  subcohort
data_for_analysis <- feather::read_feather(here::here("data/data_for_analysis.feather")) %>%
  select(omics_id, ia_subcohort)

# Join with lipidomics
temp_lips <- df_lips %>%
  left_join(data_for_analysis)


# Scale
df_lips_scaled <- temp_lips %>%
  mutate(across(-c(omics_id, ia_subcohort), ~scale_subcohort(temp_lips, .)))


# Get weights from trial (unadjusted for class sums)
mufa_lipids <- df_lips_scaled %>%
  select(all_of(mufa_weights$lipid))

pufa_lipids <- df_lips_scaled %>%
  select(all_of(pufa_weights$lipid))



# Calculate scores by multiplying weights with respective columns
diet_scores <- tibble(omics_id = df_lips$omics_id,
                      mufa_score = as.numeric(as.matrix(mufa_lipids) %*% mufa_weights$estimate),
                      pufa_score = as.numeric(as.matrix(pufa_lipids) %*% pufa_weights$estimate))


write_rds(diet_scores, here::here("data/separate_scores.rds"))

















