# pak::pkg_install("netCoupler/netCoupler")

# Read packages
library(tidyverse)
library(NetCoupler)
library(survival)
library(future)

# Read lipid data for network generation
lips <- read_rds(here::here("data/data_for_clusters.rds"))

# Get lipid names
lips_names <- lips |>
  select(-c(omics_id, ia_subcohort)) |>
  names()


# Read data for outcome analysis
data_for_analysis <- read_rds(here::here("data/data_for_analysis.rds"))

# Wrangle data (no log-transform needed, was already done!!)
pheno_df <- data_for_analysis |>
  mutate(across(c(epic_cvd, case_diab_caco, fasting, educc3, smk_4cat, alccat,
                  antihyp, lipidlower, ass), as_factor)) |>
  left_join(lips) |>
  mutate(across(all_of(lips_names), ~scale(.) |>
                  as.numeric()))


# Combine phenotype data and scores to create outcome-specific datasets
diab_df <- pheno_df |>
  filter(ia_subcohort == 1 | case_diab_caco == 1,
         start_time_diab <= stop_time_diab,
         !is.na(start_time_diab)) |>
  mutate(case_diab_caco =  if_else(case_diab_caco == 1, 1, 0)) |>
  rename(start = start_time_diab,
         stop = stop_time_diab,
         outcome = case_diab_caco)


cvd_df <- pheno_df |>
  filter(ia_subcohort == 1 | epic_cvd == 1,
         start_time_cvd <= stop_time_cvd,
         !is.na(start_time_cvd)) |>
  mutate(epic_cvd =  if_else(epic_cvd == 1, 1, 0)) |>
  rename(start = start_time_cvd,
         stop = stop_time_cvd,
         outcome = epic_cvd)

# Generate network
metabolite_network <- lips |>
  filter(ia_subcohort == 1) |>
  nc_estimate_network(all_of(lips_names))


# First iteration ---------------------------------------------------------

plan(multisession)
diab_netcoupler <- diab_df |>
  mutate(surv_object = Surv(time = start,
                            time2 = stop,
                            event = outcome)) |>
  nc_estimate_outcome_links(edge_tbl = as_edge_tbl(metabolite_network),
                            outcome = "surv_object",
                            adjustment_vars = c("strata(as.integer(age))", "sex", "height", "waist", "aktiv", "fasting", "antihyp",
                                                "lipidlower", "ass", "gj", "smk_4cat", "alccat", "educc3", "syst", "diast", "cluster(omics_id)"),
                            model_function = survival::coxph,
                            exponentiate = TRUE)



cvd_netcoupler <- cvd_df |>
  mutate(surv_object = Surv(time = start,
                            time2 = stop,
                            event = outcome)) |>
  nc_estimate_outcome_links(
    edge_tbl = as_edge_tbl(metabolite_network),
    outcome = "surv_object",
    adjustment_vars = c("strata(as.integer(age))", "sex", "height", "waist", "aktiv", "fasting", "antihyp",
                        "lipidlower", "ass", "gj", "smk_4cat", "alccat", "educc3", "syst", "diast" ,
                        "prevalent_diab", "cluster(omics_id)"),
    model_function = survival::coxph,
    exponentiate = TRUE)

plan(sequential)

effect_cvd_1 <- cvd_netcoupler |>
  filter(effect == "direct") |>
  pull(index_node)

effect_cvd_1
# "pi_fa16_0"

effect_diab_1 <- diab_netcoupler |>
  filter(effect == "direct") |>
  pull(index_node)

effect_diab_1
# "dag_fa14_0"  "dcer_fa20_0" "lcer_fa14_0" "lpc_fa14_0"


# Second iteration --------------------------------------------------------

plan(multisession)

diab_netcoupler_2 <- diab_df |>
  mutate(surv_object = Surv(time = start,
                            time2 = stop,
                            event = outcome)) |>
  nc_estimate_outcome_links(edge_tbl = as_edge_tbl(tidygraph::filter(metabolite_network, !name %in% effect_diab_1)),
                            outcome = "surv_object",
                            adjustment_vars = c("strata(as.integer(age))", "sex", "height", "waist", "aktiv", "fasting", "antihyp",
                                                "lipidlower", "ass", "gj", "smk_4cat", "alccat", "educc3", "syst", "diast",
                                                "cluster(omics_id)", effect_diab_1),
                            model_function = survival::coxph,
                            exponentiate = TRUE)



cvd_netcoupler_2 <- cvd_df |>
  mutate(surv_object = Surv(time = start,
                            time2 = stop,
                            event = outcome)) |>
  nc_estimate_outcome_links(edge_tbl = as_edge_tbl(tidygraph::filter(metabolite_network, !name %in% effect_cvd_1)),
                            outcome = "surv_object",
                            adjustment_vars = c("strata(as.integer(age))", "sex", "height", "waist", "aktiv", "fasting", "antihyp",
                                                "lipidlower", "ass", "gj", "smk_4cat", "alccat", "educc3", "syst", "diast" ,
                                                "prevalent_diab", "cluster(omics_id)", effect_cvd_1),
                            model_function = survival::coxph,
                            exponentiate = TRUE)

plan(sequential)


effect_cvd_2 <- cvd_netcoupler_2 |>
  filter(effect == "direct") |>
  pull(index_node)

effect_cvd_2
# [1] "cer_fa26_0"


effect_diab_2 <- diab_netcoupler_2 |>
  filter(effect == "direct") |>
  pull(index_node)

effect_diab_2
# [1] "ce_fa15_0"  "dag_fa14_1" "pc_fa14_0"

# Third iteration --------------------------------------------------------

plan(multisession)
diab_netcoupler_3 <- diab_df |>
  mutate(surv_object = Surv(time = start,
                            time2 = stop,
                            event = outcome)) |>
  nc_estimate_outcome_links(edge_tbl = as_edge_tbl(tidygraph::filter(metabolite_network,
                                                                     !name %in% c(effect_diab_1, effect_diab_2))),
                            outcome = "surv_object",
                            adjustment_vars = c("strata(as.integer(age))", "sex", "height", "waist", "aktiv", "fasting", "antihyp",
                                                "lipidlower", "ass", "gj", "smk_4cat", "alccat", "educc3", "syst", "diast", "cluster(omics_id)", effect_diab_1, effect_diab_2),
                            model_function = survival::coxph,
                            exponentiate = TRUE)


cvd_netcoupler_3 <- cvd_df |>
  mutate(surv_object = Surv(time = start,
                            time2 = stop,
                            event = outcome)) |>
  nc_estimate_outcome_links(edge_tbl = as_edge_tbl(tidygraph::filter(metabolite_network,
                                                                     !name %in% c(effect_cvd_1, effect_cvd_2))),
                            outcome = "surv_object",
                            adjustment_vars = c("strata(as.integer(age))", "sex", "height", "waist", "aktiv", "fasting", "antihyp",
                                                "lipidlower", "ass", "gj", "smk_4cat", "alccat", "educc3", "syst", "diast" ,
                                                "prevalent_diab", "cluster(omics_id)", effect_cvd_1, effect_cvd_2),
                            model_function = survival::coxph,
                            exponentiate = TRUE)

plan(sequential)

effect_cvd_3 <- cvd_netcoupler_3 |>
  filter(effect == "direct") |>
  pull(index_node)

effect_cvd_3
# NULL


effect_diab_3 <- diab_netcoupler_3 |>
  filter(effect == "direct") |>
  pull(index_node)

effect_diab_3
# [1] "pc_fa17_0"



plan(multisession)

diab_netcoupler_4 <- diab_df |>
  mutate(surv_object = Surv(time = start,
                            time2 = stop,
                            event = outcome)) |>
  nc_estimate_outcome_links(edge_tbl = as_edge_tbl(tidygraph::filter(metabolite_network,
                                                                     !name %in% c(effect_diab_1, effect_diab_2, effect_diab_3))),
                            outcome = "surv_object",
                            adjustment_vars = c("strata(as.integer(age))", "sex", "height", "waist", "aktiv", "fasting", "antihyp",
                                                "lipidlower", "ass", "gj", "smk_4cat", "alccat", "educc3", "syst", "diast", "cluster(omics_id)", effect_diab_1, effect_diab_2, effect_diab_3),
                            model_function = survival::coxph,
                            exponentiate = TRUE)


plan(sequential)


effect_diab_4 <- diab_netcoupler_4 |>
  filter(effect == "direct") |>
  pull(index_node)

effect_diab_4
# NULL


netcoupler_direct_effects <- list(diab = c(effect_diab_1, effect_diab_2, effect_diab_3),
                                  cvd = c(effect_cvd_1, effect_cvd_2))


bind_rows(diab_netcoupler |>
            mutate(iter = 1, outcome = "T2D"),
          diab_netcoupler_2 |>
            mutate(iter = 2, outcome = "T2D"),
          diab_netcoupler_3 |>
            mutate(iter = 3, outcome = "T2D"),
          diab_netcoupler_4 |>
            mutate(iter = 4, outcome = "T2D"),
          cvd_netcoupler |>
            mutate(iter = 0, outcome = "CVD"),
          cvd_netcoupler_2 |>
            mutate(iter = 1, outcome = "CVD"),
          cvd_netcoupler_3 |>
            mutate(iter = 2, outcome = "CVD")) |>
  xlsx::write.xlsx(file = here::here("doc/netcoupler_res.xlsx"))

