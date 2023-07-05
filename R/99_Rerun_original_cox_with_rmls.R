library(tidyverse)

data_for_analysis <- read_rds(here::here("data/data_for_analysis.rds"))

combined_score <- read_rds(here::here("data/predimed_score.rds")) %>%
  rename(predimed_score = comb_score)




# make outcome dfs
# Prep phenotype data
pheno_df <- data_for_analysis %>%
  mutate(across(c(epic_cvd, case_diab_caco, fasting, educc3, smk_4cat, alccat,
                  antihyp, lipidlower, ass), as_factor)) %>%
  left_join(combined_score) %>%
  mutate(across(c(
    "corr_trigly",
    "corr_chol",
    "corr_hdl",
    "corr_hscrp",
    "non_hdl"), ~log(.) %>%
      #scale() %>%
      as.numeric()),
    comb_score_non_hdl = (comb_score + (non_hdl * -0.11))/3.48,
    comb_score = comb_score/3.48,
    non_hdl = non_hdl / -0.11)  # this is the effect seen in DIVAS




# Combine phenotype data and scores to create outcome-specific datasets
diab_df <- pheno_df %>%
  filter(ia_subcohort == 1 | case_diab_caco == 1,
         start_time_diab <= stop_time_diab,
         !is.na(start_time_diab)) %>%
  mutate(case_diab_caco =  if_else(case_diab_caco == 1, 1, 0),
         comb_score = predimed_score) %>%
  rename(start = start_time_diab,
         stop = stop_time_diab,
         outcome = case_diab_caco)

cvd_df <- pheno_df %>%
  filter(ia_subcohort == 1 | epic_cvd == 1,
         start_time_cvd <= stop_time_cvd,
         !is.na(start_time_cvd)) %>%
  mutate(epic_cvd =  if_else(epic_cvd == 1, 1, 0),
         comb_score = predimed_score) %>%
  rename(start = start_time_cvd,
         stop = stop_time_cvd,
         outcome = epic_cvd)


# Models
# 1 age, sex
# 2 1 + height,waist, aktiv, fasting, antihyp, lipidlower, ass, gj, smk_4cat, alccat, educc3, syst, diast
# 3 wie 2 nur sex*score
# 4 2 + trigly
# 5 2 + chol
# 6 2 + hdl
# 7 2 + non-hdl
# 8 2 + hdl, chol, trigly
# 9 2 + hba1c
# 10 2 + hdl, chol, trigly, hba1c
# 10 2 + log(hscrp)



# define models as list
mods_diab_comb <- list(
  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)),

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast,

  Surv(start, stop, outcome) ~
    comb_score*sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast,

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_trigly,

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_chol,

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hdl,

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + non_hdl,

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_chol + corr_trigly + corr_hdl,

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hba1c,

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hba1c +
    corr_chol + corr_trigly + corr_hdl,

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hscrp,

  Surv(start, stop, outcome) ~
    comb_score_non_hdl + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hscrp,

  Surv(start, stop, outcome) ~
    non_hdl + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hscrp) %>%
  set_names(paste("diab_comb", 1:length(.)))


mods_cvd_comb <- list(
  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)),

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + prevalent_diab,

  Surv(start, stop, outcome) ~
    comb_score*sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + prevalent_diab,

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_trigly + prevalent_diab,

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_chol + prevalent_diab,

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hdl + prevalent_diab,

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + non_hdl + prevalent_diab,

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_chol + corr_trigly + corr_hdl +
    prevalent_diab,

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hba1c + prevalent_diab,

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hba1c +
    corr_chol + corr_trigly + corr_hdl + prevalent_diab,

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hscrp + prevalent_diab,

  Surv(start, stop, outcome) ~
    comb_score_non_hdl + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hscrp + prevalent_diab,

  Surv(start, stop, outcome) ~
    non_hdl + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + corr_hscrp + prevalent_diab) %>%
  set_names(paste("cvd_comb", 1:length(.)))




# Run models
library(survival)

diab_res <- map(mods_diab_comb, ~coxph(.x ,
                                       ties = "efron",
                                       data = diab_df))

cvd_res <- map(mods_cvd_comb, ~coxph(.x ,
                                     ties = "efron",
                                     data = cvd_df))


# collect relevant coefficients
coeffs <- c(diab_res, cvd_res) %>%
  map_dfr(~broom::tidy(.x,
                       exponentiate = TRUE,
                       conf.int = TRUE,
                       conf.level = 0.95),
          .id = "model")

# check sex interaction




temp <- coeffs %>%
  filter(str_detect(term, "score")) %>%
  separate(model, into = c("outcome", "score", "model")) %>%
  mutate(outcome = if_else(outcome == "diab", "Type 2 diabetes", "CVD"),
         model = as.double(model) ,
         model_names = case_when(
           model == 1 ~ "Age & Sex",
           model == 2 ~ "Multivariable adjusted (MV)",
           model == 3 ~ "MV interact sex",
           model == 4 ~ "MV + TG",
           model == 5 ~ "MV + TC",
           model == 6 ~ "MV + HDL-C",
           model == 7 ~ "MV + Non-HDL-C",
           model == 8 ~ "MV + HDL + TC + TG",
           model == 9 ~ "MV + HbA1c",
           model == 10 ~ "MV + HDL + TC + TG + HbA1c",
           model == 11 ~ "MV + hsCRP",
           model == 12 ~ "non_hdl_lipidomics_score",
           model == 13 ~ "Non-HDL-C")) %>%
  filter(!model %in% c(1, 3, 8, 10, 12, 13))





cox_figure <- temp %>%
  mutate(panel = if_else(model == 2, "one", "two"),
         label = paste0(sprintf('%.2f', round(estimate, 2)),
                        " (",
                        sprintf('%.2f', round(conf.low, 2)),
                        ", ",
                        sprintf('%.2f', round(conf.high, 2)),
                        ")")) %>%
  ggplot() +
  ggforestplot::geom_stripes(aes(y = fct_reorder(model_names, -model))) +
  geom_vline(xintercept = 1,
             lty = 2) +
  labs(x = NULL,
       y = NULL,
       fill = NULL) +
  theme_light() +
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
                      y = fct_reorder(model_names, -model),
                      xmin = conf.low,
                      xmax = conf.high),
                  orientation = "y",
                  shape = "square") +
  facet_wrap(~outcome) +
  scale_x_log10(breaks = c(seq(0.6, 1, 0.1),
                           map_dbl(seq(0.6, 1, 0.1), ~1/.) %>% round(digits = 1))) +
  labs(x = "HR of theoretical intervention effect (95%-CI)",
       y = NULL) +
  geom_text(inherit.aes = FALSE,
            aes(y = fct_reorder(model_names, -model),
                x = 0.3,
                label = label),
            hjust = 0)




reduction_fig <- coeffs %>%
  separate(model, into = c("outcome", "score", "model")) %>%
  mutate(outcome = if_else(outcome == "diab", "Type 2 diabetes", "CVD"),
         model = as.double(model) ,
         model_names = case_when(
           model == 1 ~ "Age & Sex",
           model == 2 ~ "Multivariable adjusted (MV)",
           model == 3 ~ "MV interact sex",
           model == 4 ~ "MV + TG",
           model == 5 ~ "MV + TC",
           model == 6 ~ "MV + HDL-C",
           model == 7 ~ "MV + Non-HDL-C",
           model == 8 ~ "MV + HDL + TC + TG",
           model == 9 ~ "MV + HbA1c",
           model == 10 ~ "MV + HDL + TC + TG + HbA1c",
           model == 11 ~ "MV + hsCRP",
           model == 12 ~ "non_hdl_lipidomics_score",
           model == 13 ~ "Non-HDL-C")) %>%
  filter(model %in% c(2, 7, 12, 13),
         str_detect(term, "score|non_hdl")) %>%
  mutate(panel = case_when(model %in% c(2, 13) ~ "Separately",
                           model == 7 ~ "Mutually adjusted",
                           model == 12 ~ "Unified score") %>%
           fct_relevel(c("Separately",
                         "Mutually adjusted",
                         "Unified score")),
         term = case_when(term == "comb_score" ~ "UFA lipidomics score",
                          term == "non_hdl" ~ "Non-HDL-C",
                          term == "comb_score_non_hdl" ~ "Unified score")) %>%
  ggplot() +
  theme_light() +
  geom_pointrange(aes(y = estimate,
                      x = panel,
                      ymin = conf.low,
                      ymax = conf.high,
                      fill = term),
                  position = position_dodge(width = 0.5),
                  shape = 22,
                  fatten = 5) +
  facet_wrap(~outcome,
             scales = "free_x") +
  scale_y_log10(breaks = c(seq(0.5, 1, 0.1)),
                labels =  paste0(c(seq(0.5, 1, 0.1) - 1) * 100, "%")) +
  geom_hline(yintercept = 1,
             lty = 2) +
  labs(x = NULL,
       y = "% risk reduction\nper DIVAS intervention effect",
       fill = NULL) +
  theme_light() +
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
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(vjust = 0)) +
  ggsci::scale_fill_futurama()



library(ggpubr)


figure <- ggarrange(cox_figure, reduction_fig, ncol = 1,
                    labels = "AUTO",
                    heights = c(1, 0.7))
figure

# Set path
path <- here::here("doc", "Figures", "figure_4_predimed")

# Save as pdf
ggsave(plot = figure,
       glue::glue("{path}.pdf"),
       device = cairo_pdf,
       width = 22,
       height = 18,
       units = "cm")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)


