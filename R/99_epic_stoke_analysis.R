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
    mls = comb_score) # this is the effect seen in DIVAS


rmls_df <- read_rds(here::here("data/predimed_score.rds")) %>%
  rename(rmls = comb_score)


# Combine phenotype data and scores to create outcome-specific datasets
stroke_df <- pheno_df %>%
  filter(ia_subcohort == 1 | epic_stroke == 1,
         start_time_stroke <= stop_time_stroke,
         !is.na(start_time_stroke)) %>%
  mutate(epic_stroke = if_else(epic_stroke == 1, 1, 0)) %>%
  rename(start = start_time_stroke,
         stop = stop_time_stroke,
         outcome = epic_stroke) %>%
  left_join(rmls_df) %>%
  mutate(across(c(mls, rmls), ~log(.) %>% scale() %>% as.numeric()))



# define models as list
mods_stroke_mls <- list(
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
    smk_4cat + alccat + educc3 + syst + diast + corr_hscrp + prevalent_diab)


mods_stroke_rmls <- list(
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
    smk_4cat + alccat + educc3 + syst + diast + corr_hscrp + prevalent_diab)




# Run models
library(survival)




stroke_mls_res <- map(mods_stroke_mls, ~coxph(.x ,
                                     ties = "efron",
                                     data = stroke_df)) %>%
  map_dfr(~broom::tidy(.x,
                       exponentiate = TRUE,
                       conf.int = TRUE,
                       conf.level = 0.95),
          .id = "model") %>%
  filter(term == "mls") %>%
  mutate(sensi_w = "no")

stroke_mls_res_w <- map(mods_stroke_mls, ~coxph(.x ,
                                             ties = "efron",
                                             data = stroke_df %>%
                                               filter(sex == 2))) %>%
  map_dfr(~broom::tidy(.x,
                       exponentiate = TRUE,
                       conf.int = TRUE,
                       conf.level = 0.95),
          .id = "model") %>%
  filter(term == "mls") %>%
  mutate(sensi_w = "yes")


stroke_rmls_res <- map(mods_stroke_rmls, ~coxph(.x ,
                                              ties = "efron",
                                              data = stroke_df)) %>%
  map_dfr(~broom::tidy(.x,
                       exponentiate = TRUE,
                       conf.int = TRUE,
                       conf.level = 0.95),
          .id = "model") %>%
  filter(term == "rmls") %>%
  mutate(sensi_w = "no")

stroke_rmls_res_w <- map(mods_stroke_rmls, ~coxph(.x ,
                                                ties = "efron",
                                                data = stroke_df %>%
                                                  filter(sex == 2))) %>%
  map_dfr(~broom::tidy(.x,
                       exponentiate = TRUE,
                       conf.int = TRUE,
                       conf.level = 0.95),
          .id = "model") %>%
  filter(term == "rmls") %>%
  mutate(sensi_w = "yes")






temp <- bind_rows(stroke_mls_res, stroke_mls_res_w, stroke_rmls_res, stroke_rmls_res_w) %>%
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
         score = case_when(term == "mls" & sensi_w == "no" ~ "MLS all",
                           term == "mls" & sensi_w == "yes" ~ "MLS only women",
                           term == "rmls" & sensi_w == "no" ~ "rMLS all",
                           term == "rmls" & sensi_w == "yes" ~ "rMLS only women"))

# Save temporary results
temp %>%
  select(-model) %>%
  select(term, model_names, sensi_w, everything()) %>%
  write_tsv(here::here("doc/Tables/cox_stroke_mls_rmls.tsv"))

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

figure

# Set path
path <- here::here("doc", "Figures", "stroke_mls_rmls_women")

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


