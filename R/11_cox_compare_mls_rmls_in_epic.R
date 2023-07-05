library(tidyverse)

data_for_analysis <- read_rds(here::here("data/data_for_analysis.rds"))

# make outcome dfs
# Prep phenotype data
pheno_df <- data_for_analysis %>%
  mutate(across(c(epic_cvd, case_diab_caco, fasting, educc3, smk_4cat, alccat,
                  antihyp, lipidlower, ass), as_factor)) %>%
  mutate(across(c(
    "corr_trigly",
    "corr_chol",
    "corr_hdl",
    "corr_hscrp",
    "non_hdl"), ~log(.) %>%
      #scale() %>%
      as.numeric()),
    mls = comb_score)


rmls_df <- read_rds(here::here("data/predimed_score.rds")) %>%
  rename(rmls = comb_score)


# Combine phenotype data and scores to create outcome-specific datasets
diab_df <- pheno_df %>%
  filter(ia_subcohort == 1 | case_diab_caco == 1,
         start_time_diab <= stop_time_diab,
         !is.na(start_time_diab)) %>%
  mutate(case_diab_caco =  if_else(case_diab_caco == 1, 1, 0)) %>%
  rename(start = start_time_diab,
         stop = stop_time_diab,
         outcome = case_diab_caco) %>%
  left_join(rmls_df) %>%
  mutate(across(c(mls, rmls), ~log(.) %>% scale() %>% as.numeric()))

cvd_df <- pheno_df %>%
  filter(ia_subcohort == 1 | epic_cvd == 1,
         start_time_cvd <= stop_time_cvd,
         !is.na(start_time_cvd)) %>%
  mutate(epic_cvd =  if_else(epic_cvd == 1, 1, 0)) %>%
  rename(start = start_time_cvd,
         stop = stop_time_cvd,
         outcome = epic_cvd) %>%
  left_join(rmls_df) %>%
  mutate(across(c(mls, rmls), ~log(.) %>% scale() %>% as.numeric()))




# define models as list
mods_diab_mls <- list(
  Surv(start, stop, outcome) ~
    mls + sex + cluster(omics_id) + strata(round(age, 0)),

  Surv(start, stop, outcome) ~
    mls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast,

  Surv(start, stop, outcome) ~
    mls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_trigly,

  Surv(start, stop, outcome) ~
    mls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_chol,

  Surv(start, stop, outcome) ~
    mls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hdl,

  Surv(start, stop, outcome) ~
    mls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + non_hdl,

  Surv(start, stop, outcome) ~
    mls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hba1c,

  Surv(start, stop, outcome) ~
    mls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hscrp) %>%
  set_names(paste("diab_mls", 1:length(.)))


mods_cvd_mls <- list(
  Surv(start, stop, outcome) ~
    mls + sex + cluster(omics_id) + strata(round(age, 0)),

  Surv(start, stop, outcome) ~
    mls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + prevalent_diab,

  Surv(start, stop, outcome) ~
    mls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_trigly + prevalent_diab,

  Surv(start, stop, outcome) ~
    mls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_chol + prevalent_diab,

  Surv(start, stop, outcome) ~
    mls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hdl + prevalent_diab,

  Surv(start, stop, outcome) ~
    mls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + non_hdl + prevalent_diab,

  Surv(start, stop, outcome) ~
    mls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hba1c + prevalent_diab,

  Surv(start, stop, outcome) ~
    mls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hscrp + prevalent_diab) %>%
  set_names(paste("cvd_mls", 1:length(.)))


mods_diab_rmls <- list(
  Surv(start, stop, outcome) ~
    rmls + sex + cluster(omics_id) + strata(round(age, 0)),

  Surv(start, stop, outcome) ~
    rmls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast,

  Surv(start, stop, outcome) ~
    rmls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_trigly,

  Surv(start, stop, outcome) ~
    rmls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_chol,

  Surv(start, stop, outcome) ~
    rmls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hdl,

  Surv(start, stop, outcome) ~
    rmls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + non_hdl,

  Surv(start, stop, outcome) ~
    rmls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hba1c,

  Surv(start, stop, outcome) ~
    rmls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hscrp) %>%
  set_names(paste("diab_rmls", 1:length(.)))


mods_cvd_rmls <- list(
  Surv(start, stop, outcome) ~
    rmls + sex + cluster(omics_id) + strata(round(age, 0)),

  Surv(start, stop, outcome) ~
    rmls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + prevalent_diab,

  Surv(start, stop, outcome) ~
    rmls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_trigly + prevalent_diab,

  Surv(start, stop, outcome) ~
    rmls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_chol + prevalent_diab,

  Surv(start, stop, outcome) ~
    rmls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hdl + prevalent_diab,

  Surv(start, stop, outcome) ~
    rmls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + non_hdl + prevalent_diab,

  Surv(start, stop, outcome) ~
    rmls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hba1c + prevalent_diab,

  Surv(start, stop, outcome) ~
    rmls + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hscrp + prevalent_diab) %>%
  set_names(paste("cvd_rmls", 1:length(.)))



# Run models
library(survival)

diab_res <- map(c(mods_diab_mls, mods_diab_rmls), ~coxph(.x ,
                                       ties = "efron",
                                       data = diab_df))

cvd_res <- map(c(mods_cvd_mls, mods_cvd_rmls), ~coxph(.x ,
                                     ties = "efron",
                                     data = cvd_df))



temp <- c(cvd_res, diab_res) %>%
  map_dfr(~broom::tidy(.x,
                       exponentiate = TRUE,
                       conf.int = TRUE,
                       conf.level = 0.95),
          .id = "model") %>%
  filter(str_detect(term, "mls$")) %>%
  separate(model, into = c("outcome", "score", "model")) %>%
  mutate(model = as.double(model) ,
         model_names = case_when(
           model == 1 ~ "Age & Sex",
           model == 2 ~ "Multivariable adjusted (MV)",
           model == 3 ~ "MV + TG",
           model == 4 ~ "MV + TC",
           model == 5 ~ "MV + HDL-C",
           model == 6 ~ "MV + Non-HDL-C",
           model == 7 ~ "MV + HbA1c",
           model == 8 ~ "MV + hsCRP"),
         score = str_c(str_to_upper(outcome), " ", str_to_upper(score)) %>%
           str_replace("RMLS", "rMLS") %>%
           str_replace("DIAB", "T2D"))

# Save temporary results
temp %>%
  select(-c(model, term)) %>%
  select(outcome, score, model_names, everything()) %>%
  write_tsv(here::here("doc/Tables/cox_mls_rmls.tsv"))

figure <- temp %>%
  mutate(label = paste0(sprintf('%.2f', round(estimate, 2)),
                        " (",
                        sprintf('%.2f', round(conf.low, 2)),
                        ", ",
                        sprintf('%.2f', round(conf.high, 2)),
                        ")")) %>%
  ggplot() +
  ggforestplot::geom_stripes(aes(y = fct_rev(score))) +
  geom_vline(xintercept = 1,
             lty = 2) +
  labs(x = NULL,
       y = NULL,
       fill = NULL) +
  theme_light(base_family = "RobotoCondensed-Regular") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1),
        panel.grid = element_line(color = "grey70"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        axis.ticks.y = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = "#333333"),
        axis.text.y = element_text(hjust = 0),
        strip.background = element_rect(fill = "#333333"),
        panel.grid.major.y = element_blank()) +
  geom_pointrange(aes(x = estimate,
                      y = fct_rev(score),
                      xmin = conf.low,
                      xmax = conf.high),
                  orientation = "y",
                  shape = "square") +
  facet_wrap(~fct_reorder(model_names, model),
             ncol = 1,
             scales = "free_y") +
  scale_x_log10(breaks = c(c(0.6, 0.8, 1),
                           map_dbl(c(0.6, 0.8, 1), ~1/.) %>% round(digits = 1))) +
  labs(x = "HR per SD (95%-CI)",
       y = NULL) +
  geom_text(inherit.aes = FALSE,
            aes(y = fct_rev(score),
                x = 0.3,
                label = label),
            hjust = 0)



# Set path
path <- here::here("doc", "Figures", "mls_vs_rmls")

# Save as pdf
ggsave(plot = figure,
       glue::glue("{path}.pdf"),
       device = cairo_pdf,
       width = 22,
       height = 20,
       units = "cm")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)


