library(tidyverse)

data_for_analysis <- feather::read_feather(here::here("data/data_for_analysis.feather"))

# Prep phenotype data
pheno_df <- data_for_analysis %>%
  mutate(across(c(epic_cvd, case_diab_caco, fasting, educc3, smk_4cat, alccat,
                  antihyp, lipidlower, ass), as_factor))




# Combine phenotype data and scores to create outcome-specific datasets
diab_df <- pheno_df %>%
  filter(ia_subcohort == 1 | case_diab_caco == 1,
         start_time_diab <= stop_time_diab,
         !is.na(start_time_diab)) %>%
  mutate(case_diab_caco =  if_else(case_diab_caco == 1, 1, 0)) %>%
  rename(start = start_time_diab,
         stop = stop_time_diab,
         outcome = case_diab_caco)

cvd_df <- pheno_df %>%
  filter(ia_subcohort == 1 | epic_cvd == 1,
         start_time_cvd <= stop_time_cvd,
         !is.na(start_time_cvd)) %>%
  mutate(epic_cvd =  if_else(epic_cvd == 1, 1, 0)) %>%
  rename(start = start_time_cvd,
         stop = stop_time_cvd,
         outcome = epic_cvd)



# define models as list
mods_diab <- list(

  Surv(start, stop, outcome) ~
    cluster(omics_id) + strata(round(age, 0)) +
    sex + height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + mufa_score_divas,

  Surv(start, stop, outcome) ~
    cluster(omics_id) + strata(round(age, 0)) +
    sex + height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + non_hdl_mufa,

  Surv(start, stop, outcome) ~
    cluster(omics_id) + strata(round(age, 0)) +
    sex + height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + non_hdl_mufa + mufa_score_divas,

  Surv(start, stop, outcome) ~
    cluster(omics_id) + strata(round(age, 0)) +
    sex + height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + pufa_score_divas,

  Surv(start, stop, outcome) ~
    cluster(omics_id) + strata(round(age, 0)) +
    sex + height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + non_hdl_pufa,

  Surv(start, stop, outcome) ~
    cluster(omics_id) + strata(round(age, 0)) +
    sex + height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + non_hdl_pufa + pufa_score_divas) %>%
  set_names(paste("diab", 1:length(.)))


mods_cvd <- list(

  Surv(start, stop, outcome) ~
    cluster(omics_id) + strata(round(age, 0)) +
    sex + height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + prevalent_diab + mufa_score_divas,

  Surv(start, stop, outcome) ~
    cluster(omics_id) + strata(round(age, 0)) +
    sex + height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + prevalent_diab + non_hdl_mufa,

  Surv(start, stop, outcome) ~
    cluster(omics_id) + strata(round(age, 0)) +
    sex + height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + prevalent_diab + non_hdl_mufa + mufa_score_divas,


  Surv(start, stop, outcome) ~
    cluster(omics_id) + strata(round(age, 0)) +
    sex + height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + prevalent_diab + pufa_score_divas,

  Surv(start, stop, outcome) ~
    cluster(omics_id) + strata(round(age, 0)) +
    sex + height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + prevalent_diab + non_hdl_pufa,

  Surv(start, stop, outcome) ~
    cluster(omics_id) + strata(round(age, 0)) +
    sex + height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + prevalent_diab + non_hdl_pufa + pufa_score_divas) %>%
  set_names(paste("cvd", 1:length(.)))









# Run models
library(survival)

diab_res <- map(c(mods_diab) , ~coxph(.x ,
                                      ties = "efron",
                                      data = diab_df))

cvd_res <- map(c(mods_cvd), ~coxph(.x ,
                                   ties = "efron",
                                   data = cvd_df))



# collect relevant coefficients
coeffs <- c(diab_res, cvd_res) %>%
  map_dfr(~broom::tidy(.x,
                       exponentiate = TRUE,
                       conf.int = TRUE,
                       conf.level = 0.95),
          .id = "model") %>%
  filter(term %in% c("non_hdl_mufa", "non_hdl_pufa", "mufa_score_divas", "pufa_score_divas")) %>%
  separate(model, into = c("outcome", "model")) %>%
  select(outcome, term, model, estimate,   p.value, conf.low, conf.high, std.error, robust.se, statistic) %>%
  mutate(term = str_remove(term, "corr_"),
         diet = if_else(str_detect(term, "mufa"), "MUFA", "Mixed UFA") %>%
           fct_rev(),
         model = if_else(model %in% c(3, 6), "Mutually adjusted", "Individually") %>% fct_rev(),
         term = if_else(str_detect(term, "non_hdl"), "Non-HDL-C", "Lipidomics Score"),
         outcome = if_else(outcome == "diab", "Type 2 diabetes", "CVD")) %>%
  select(outcome, diet, term, model, everything()) %>%
  arrange(outcome, diet, model)



write_tsv(coeffs, here::here("doc/Tables/compare_non_hdl.tsv"))


figure <- coeffs %>%
  mutate(model = str_replace(model, "Mutually adjusted", "Mutually\nadjusted") %>%
           fct_rev()) %>%
  ggplot() +
  ggforestplot::geom_stripes(aes(y = model)) +
  geom_pointrange(aes(y = model,
                      x = estimate,
                      xmin = conf.low,
                      xmax = conf.high,
                      fill = term),
                  shape = 21,
                  position = position_dodge2(width = 0.5)) +
  facet_grid(diet~outcome,
             scales = "free_y") +
  geom_vline(xintercept = 1,
             lty = 2) +
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
  coord_cartesian(clip = "off") +
  scale_x_log10() +
  scale_fill_manual(values = c("Lipidomics Score" = "#ff0000",
                               "Non-HDL-C" = "#4a62ff")) +
  labs(x = "HR of theoretical intervention effect (95%-CI)",
       y = NULL,
       fill = NULL)


# Set path
path <- here::here("doc", "Figures", "compare_non_hdl")

# Save as pdf
ggsave(plot = figure,
       glue::glue("{path}.pdf"),
       device = cairo_pdf,
       width = 20,
       height = 10,
       units = "cm")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)


