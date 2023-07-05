library(tidyverse)
source(here::here("R/00_functions.R"))


# Read DIVAS data
divas <- read_rds(here::here("data/Within_class_FA_PE.rds")) %>%
  janitor::clean_names() %>%
  mutate(diet = if_else(diet == "A", "SFA", "UFA"),
         across(contains("_v"), ~log(.) %>%
                  #scale() %>%
                  as.numeric())) %>%
  rename(sex = gender)


divas_z <- read_rds(here::here("data/Within_class_FA_PE.rds")) %>%
  janitor::clean_names() %>%
  mutate(diet = if_else(diet == "A", "SFA", "UFA"),
         across(contains("_v"), ~log(.) %>%
                  scale() %>%
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


matched_epic_predimed <- read_rds(here::here("matched_epic_predimed.rds")) %>%
  filter(!is.na(metabolite)) %>%
  select(within_class_fa_sum) %>%
  unique()


weights_matched <- weights %>%
  filter(lipid_display %in% matched_epic_predimed$within_class_fa_sum)

# write for clemens
weights_matched %>%
  write_tsv(here::here("data/weights_for_reduced_score.tsv"))


# Read and wrangle EPIC-Potsdam data and scale
df_lips <- read_rds(here::here("data/lips_for_predimed.rds"))


# Get weights from trial (unadjusted for class sums)
lipids_to_score <- df_lips %>%
  select(all_of(weights_matched$lipid))


# Calculate score by multiplying weights with respective columns
combined_score <- tibble(omics_id = df_lips$omics_id,
                         comb_score = as.numeric(as.matrix(lipids_to_score) %*% weights_matched$estimate))


write_rds(combined_score, here::here("data/predimed_score.rds"))



