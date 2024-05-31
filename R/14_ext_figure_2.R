# Functions and packages
library(tidyverse)
source(here::here("R/00_functions.R"))


# Read DIVAS data on clinical markers
divas_trad <- readxl::read_xlsx(here::here("data/Lipids_DIVAS_for_Fabian.xlsx")) %>%
  janitor::clean_names() %>%
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
                                "homa_ir")) %>%
  mutate(model = map(var, ~lm.wrapper(divas_trad, .)),
         coeffs = map(model, ~broom::tidy(., conf.int = TRUE, conf.level = 0.95))) %>%
  unnest(coeffs) %>%
  filter(str_detect(term, "diet")) %>%
  mutate(p.value_adj = p.adjust(p.value, method = "fdr"))


# Pull out effect on non-hdl
weights <- divas_effects %>%
  filter(var %in% c("non_hdl", "hdl", "tag", "glucose", "crp")) %>%
  select(estimate) %>%
  deframe()


ids <- divas_trad |>
  select(volunteer, non_hdl_v1, hdl_v1, tag_v1, glucose_v1, crp_v1) |>
  na.omit()


clinic_v1 <- divas_trad |>
  select(non_hdl_v1, hdl_v1, tag_v1, glucose_v1, crp_v1) |>
  na.omit()

# for post intervention
clinic_v2 <- divas_trad |>
  select(non_hdl_v2, hdl_v2, tag_v2, glucose_v2, crp_v2) |>
  na.omit()


combined_score_divas <- read_rds(here::here("data/combined_score_divas.rds"))

all_score_divas <- divas_trad |>
  select(volunteer, diet, age, sex, bmi) |>
  filter(volunteer %in% ids$volunteer) |>
  bind_cols(tibble(clinical_score_v1 = -as.numeric(as.matrix(clinic_v1) %*% weights),
                   clinical_score_v2 = -as.numeric(as.matrix(clinic_v2) %*% weights)))


all_score_divas |>
  select(contains("score")) |>
  pivot_longer(everything()) |>
  ggplot(aes(value)) +
  geom_density() +
  facet_wrap(~name, scales = "free")

tibble(var = c("clinical_score")) %>%
  mutate(model = map(var, ~lm.wrapper(na.omit(all_score_divas), .)),
         coeffs = map(model, ~broom::tidy(., conf.int = TRUE, conf.level = 0.95))) %>%
  unnest(coeffs) %>%
  filter(str_detect(term, "diet"))


# ==> 0.0106

data_for_analysis <- read_rds(here::here("data/data_for_analysis.rds"))
ids <- data_for_analysis |>
  select(omics_id, non_hdl, corr_hdl, corr_trigly, corr_glucose, corr_hscrp) |>
  na.omit()



prepped_biomarkers <- data_for_analysis |>
  select(non_hdl, corr_hdl, corr_trigly, corr_glucose, corr_hscrp) |>
  mutate(non_hdl = convert_cholesterol(non_hdl),
         hdl = convert_cholesterol(corr_hdl),
         tag = convert_triglycerides(corr_trigly),
         glucose = convert_glucose(corr_glucose),
         crp = corr_hscrp,
         .keep = "none") |>
  mutate(across(everything(), log)) |>
  na.omit()


ready_df <- data_for_analysis |>
  filter(omics_id %in% ids$omics_id) |>
  bind_cols(tibble(clinical_score = as.numeric(as.matrix(prepped_biomarkers) %*% weights)/0.0106)) |>
  mutate(comb_score =  comb_score/3.48)


hist(ready_df$clinical_score)
hist(ready_df$comb_score)



# Combine phenotype data and scores to create outcome-specific datasets
diab_df <- ready_df %>%
  filter(ia_subcohort == 1 | case_diab_caco == 1,
         start_time_diab <= stop_time_diab,
         !is.na(start_time_diab)) %>%
  mutate(case_diab_caco =  if_else(case_diab_caco == 1, 1, 0)) %>%
  rename(start = start_time_diab,
         stop = stop_time_diab,
         outcome = case_diab_caco)

cvd_df <- ready_df %>%
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
    clinical_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast,

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast) %>%
  set_names(paste("diab", 1:length(.)))



mods_cvd <- list(

  Surv(start, stop, outcome) ~
    clinical_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + prevalent_diab,

  Surv(start, stop, outcome) ~
    comb_score + sex + cluster(omics_id) + strata(round(age, 0)) +
    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
    smk_4cat + alccat + educc3 + syst + diast + prevalent_diab) %>%
  set_names(paste("cvd", 1:length(.)))




# Run models
library(survival)

diab_res <- map(mods_diab, ~coxph(.x ,
                                  ties = "efron",
                                  data = diab_df))

cvd_res <- map(mods_cvd, ~coxph(.x ,
                                ties = "efron",
                                data = cvd_df))






# collect relevant coefficients
coeffs <- c(diab_res, cvd_res) %>%
  map_dfr(~broom::tidy(.x,
                       exponentiate = TRUE,
                       conf.int = TRUE,
                       conf.level = 0.95),
          .id = "model")



figure <- coeffs %>%
  filter(str_detect(term, "score")) |>
  separate(model, into = c("outcome", "model")) %>%
  mutate(outcome = if_else(outcome == "diab", "Type 2 diabetes", "CVD"),
         term = if_else(term == "comb_score", "MLS", "Clinical score")) |>
  ggplot() +
  geom_hline(yintercept = 1,
             lty = 2) +
  geom_pointrange(aes(y = estimate,
                      x = outcome,
                      ymin = conf.low,
                      ymax = conf.high,
                      fill = term),
                  position = position_dodge(width = 0.2),
                  shape = 21,
                  size = 0.5) +
  scale_y_log10(breaks = c(seq(0.5, 1, 0.1)),
                labels =  paste0(c(seq(0.5, 1, 0.1) - 1) * 100, "%")) +
  labs(x = NULL,
       y = "% risk reduction\nper DIVAS intervention effect\n",
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
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(vjust = 0)) +
  ggsci::scale_fill_d3()

# Set path
path_to <- ("/home/fabian/alle_shortcut/!MEP/Projekte/EPIC-Potsdam/Diabetes/Lipidomics/MultiLipidScore/AIP/")

path <- str_c(path_to, "Extended_Figure_2")
# Save as pdf
ggsave(plot = figure,
       glue::glue("{path}.pdf"),
       device = cairo_pdf,
       height = 10,
       width = 10,
       units = "cm")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)


