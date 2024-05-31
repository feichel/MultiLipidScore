library(tidyverse)


rmls <- read_tsv(here::here("doc/Tables/cox_mls_rmls.tsv")) |>
  filter(str_detect(score, "rMLS"),
         model != 1) |>
  select(outcome, model_names, estimate, conf.low, conf.high, statistic, p.value, model) |>
  mutate(outcome = if_else(outcome == "cvd", "CVD", "Type 2 diabetes"),
         score = "rMLS")


rmls
