library(tidyverse)

# get predimed variable names and harmonize with epic names
predimed_lipids <- read_csv(here::here("data/CVD_C8_pos_update_Met_annotation.csv"),
                            col_select = 2,
                            n_max = 203,
                            show_col_types = FALSE,
                            name_repair =  janitor::make_clean_names) %>%
  filter(!metabolite %in% c("C32:2 DAG", "C34:1 PI", "C36:4 PI")) %>%
  mutate(lipid = str_replace(metabolite, " plasmalogen", "P") %>%
           str_replace("Ceramide \\(d18:1\\)", "Cer") %>%
           str_replace("AG", "G") %>%
           str_remove("^C") %>%
           str_remove("-A") %>%
           str_remove("-B")) %>%
  unique() %>%
  filter(!str_detect(lipid, "hydrox|sphingo")) %>%
  separate(lipid,
           into = c("fa", "base"),
           remove = FALSE,
           sep = " ") %>%
  separate(fa, into = c("l", "db")) %>%
  mutate(lipid_display = if_else(!base %in% c("CE", "Cer", "LPC", "LPE", "MG", "SM"),
                                 str_c(base, "(", l, "_", db, ")"),
                                 str_c(base, "(", l, ":", db, ")"))) %>%
  select(metabolite, lipid_display)


# get epic names
fa_species <- read_csv("/home/fabian/lipidomics.cardiometabolic.diseases/doc/Tables/full_results.csv",
                       show_col_types = FALSE) %>%
  filter(type == "molecular species",
         !str_detect(base, "FFA")) %>%
  select(lipid) %>%
  unique()

# produce matching tables to link species to within-class Fa abundances
matching_table <- fa_species %>%
  mutate(lipid_new = str_remove_all(lipid, "\\)") %>%
           str_replace("\\(P-", "P(") %>%
           str_replace("\\(O-", "O(") %>%
           str_remove("FA")) %>%
  separate(lipid_new, into = c("base", "l1", "db1", "l2", "db2"),
           fill = "right") %>%
  filter(!is.na(l2)) %>%
  mutate(species = if_else(base != "TG",
                           str_c(base, "(", as.numeric(l1) + as.numeric(l2), "_",
                                 as.numeric(db1) + as.numeric(db2), ")"),
                           str_c(base, "(", l1, "_", db1, ")")),
         fa1 = if_else(base != "TG",
                       str_c(base, "(FA", l1, ":", db1, ")"),
                       str_c(base, "(FA", l2, ":", db2, ")")),
         fa2 = if_else(base != "TG",
                       str_c(base, "(FA", l2, ":", db2, ")"),
                       as.character(NA))) %>%
  pivot_longer(c(fa1, fa2),
               values_to = "within_fa_abundance",
               names_to = "fa")



# Load weights that were part of UFA score
comb_weights_z <- read_rds(here::here("data/weigths_comb_score.rds")) %>%
  select(lipid, lipid_display, base, l = length, db)

weights_with_matching_info <- comb_weights_z %>%
  select(lipid_display) %>%
  left_join(matching_table, by = c("lipid_display" = "within_fa_abundance"))


one_fa_lipids <- weights_with_matching_info %>%
  filter(is.na(species)) %>%
  select(lipid_display)

more_fa_lipids <- weights_with_matching_info %>%
  filter(!is.na(species)) %>%
  select(lipid_display = species, within_class_fa_sum = lipid_display)


matched_epic_predimed <- bind_rows(one_fa_lipids, more_fa_lipids) %>%
  left_join(predimed_lipids) %>%
  mutate(within_class_fa_sum = if_else(is.na(within_class_fa_sum), lipid_display, within_class_fa_sum))


# For how many original lipids do we have matches?
matched_epic_predimed %>%
  filter(!is.na(metabolite)) %>%
  count(within_class_fa_sum)


write_rds(matched_epic_predimed,
          here::here("doc/matching_for_predimed.rds"))


matched_epic_predimed %>% na.omit %>%
  group_by(within_class_fa_sum) %>%
  mutate(metabolite) %>%
  summarize(predimed = str_c(metabolite, collapse = ";"))



# get cluster info
cluster_info <- read_rds(here::here("data/cluster.info.rds")) %>%
  select(lipid, lipid_display, louvain) %>%
  unique()


# how many matches do we have per cluster
matched_epic_predimed %>%
  left_join(cluster_info, by = c("within_class_fa_sum" = "lipid_display")) %>%
  filter(!is.na(metabolite)) %>%
  select(within_class_fa_sum, louvain) %>%
  unique() %>%
  count(louvain)


matched_epic_predimed %>%
  write_rds(here::here("matched_epic_predimed.rds"))


matched_epic_predimed %>%
  xlsx::write.xlsx(., here::here("data/matched.epic.predimed.xlsx"))














