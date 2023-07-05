# Functions and packages
library(tidyverse)
source(here::here("R/00_functions.R"))

# Read scores
combined_score_divas <- read_rds(here::here("data/combined_score_divas.rds")) %>%
  select(volunteer, comb_score_v1, comb_score_v2)

# Read DIVAS data on clinical markers
divas_trad <- readxl::read_xlsx(here::here("data/Lipids_DIVAS_for_Fabian.xlsx")) %>%
  janitor::clean_names() %>%
  left_join(combined_score_divas) %>%
  mutate(across(c(total_cholesterol_v1:il6_v2), ~na_if(., 0)),
         across(c(total_cholesterol_v1:il6_v2), ~log(.) %>%
                  #scale() %>%
                  as.numeric()),
         # across(contains("comb_score"), ~scale(.) %>%
         #          as.numeric()),
         diet = if_else(diet == "A", "SFA", "UFA"))

# Calculate effects on score and clinical markers
divas_effects <- tibble(var = c("hdl",
                                "non_hdl",
                                "tag",
                                "glucose",
                                "crp",
                                "il6",
                                "homa_ir",
                                "comb_score")) %>%
  mutate(model = map(var, ~lm.wrapper(divas_trad, .)),
         coeffs = map(model, ~broom::tidy(., conf.int = TRUE, conf.level = 0.95))) %>%
  unnest(coeffs) %>%
  filter(str_detect(term, "diet")) %>%
  mutate(p.value_adj = p.adjust(p.value, method = "fdr"))


# Pull out effect on non-hdl
weight_non_hdl <- divas_effects %>%
  filter(var == "non_hdl") %>%
  select(estimate) %>%
  deframe()

# Add weighted non-hdl to MLS
divas_trad_updated <- divas_trad %>%
  mutate(comb_score_non_hdl_v1 = comb_score_v1 + (non_hdl_v1 * weight_non_hdl),
         comb_score_non_hdl_v2 = comb_score_v2 + (non_hdl_v2 * weight_non_hdl))

# Check effect on MLS-non-hdl score
divas_effects_lipidscore_non_hdl <- tibble(var = "comb_score_non_hdl") %>%
  mutate(model = map(var, ~lm.wrapper(divas_trad_updated, .)),
         coeffs = map(model, ~broom::tidy(., conf.int = TRUE, conf.level = 0.95))) %>%
  unnest(coeffs) %>%
  filter(str_detect(term, "diet")) %>%
  mutate(p.value_adj = p.adjust(p.value, method = "fdr"))

# save effects
bind_rows(divas_effects,
          divas_effects_lipidscore_non_hdl) %>%
  mutate(var = case_when(var == "hdl" ~ "HDL-C",
                         var == "non_hdl" ~ "Non-HDL-C",
                         var == "glucose" ~ "Fasting glucose",
                         var == "crp" ~ "hsCRP",
                         var == "il6" ~ "IL-6",
                         var == "homa_ir" ~ "HOMA-IR",
                         var == "comb_score" ~ "Lipidomics-score",
                         var == "comb_score_non_hdl" ~ "Lipidomics-score + non-HDL",
                         TRUE ~ "TG")) %>%
  mutate(var = fct_relevel(var, c("HDL-C", "Non-HDL-C", "TG", "Fasting glucose", "HOMA-IR",
                                  "hsCRP", "IL-6", "Lipidomics-score", "Lipidomics-score + non-HDL")))















