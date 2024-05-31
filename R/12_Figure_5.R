library(tidyverse)


# NHS data ----------------------------------------------------------------

# input tables provided

change_t2d <- read_csv(here::here("data/fabian/teny_change_MLS_and_t2D.csv")) |>
  slice(1) |>
  mutate(term = "T2D (10y-change)",
         bmi_adjusted = "BMI-adjusted")


baseline_t2d_adj <- read_csv(here::here("data/fabian/baseline_MLS_and_t2D_CI.csv")) |>
  slice(1) |>
  mutate(term = "T2D (baseline)",
         bmi_adjusted = "BMI-adjusted")

baseline_t2d <- read_csv(here::here("data/fabian/baseline_MLS_and_t2D_CI_unadj.csv")) |>
  mutate(term = "T2D (baseline)",
         bmi_adjusted = "BMI-unadjusted")


baseline_stroke <- read_csv(here::here("data/fabian/baseline_MLS_and_stroke.csv")) |>
  mutate(term = "Stroke (baseline)",
         bmi_adjusted = "BMI-unadjusted")
baseline_stroke_adj <- read_csv(here::here("data/fabian/baseline_MLS_and_stroke_adj.csv")) |>
  slice(1) |>
  mutate(term = "Stroke (baseline)",
         bmi_adjusted = "BMI-adjusted")




nhs <- bind_rows(
  change_t2d,
  baseline_t2d_adj,
  baseline_t2d,
  baseline_stroke,
  baseline_stroke_adj) |>
  janitor::clean_names() |>
  rename(outcome = term) |>
  ggplot() +
  geom_vline(xintercept = 1) +
  geom_pointrange(aes(y = fct_rev(fct_relevel(outcome,
                                              c("Stroke (baseline)",
                                                "T2D (baseline)",
                                                "T2D (10y-change)"))),
                      x = hr,
                      xmin = lcl,
                      xmax = ucl,
                      fill = bmi_adjusted),
                  shape = 22,
                  fatten = 8,
                  position = position_dodge(width = 0.5)) +
  ggforestplot::geom_stripes(aes(y = outcome),
                             show.legend = FALSE) +
  theme_light(base_family = "RobotoCondensed-Regular") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major = element_line(size = 0.1),
        panel.grid = element_line(color = "grey70"),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(color = "#333333"),
        panel.grid.major.y = element_blank(),
        legend.position = "bottom",
        axis.ticks = element_blank()) +
  labs(y = NULL,
       fill = NULL,
       x = "Odds ratio per SD increase in rMLS (95%-CI)") +
  ggsci::scale_fill_npg(breaks = c("BMI-unadjusted", "BMI-adjusted")) +
  scale_x_log10(breaks = c(0.6, 0.7,0.8, 0.9, 1), limits = c(0.4, 1.05)) +
  geom_text(inherit.aes = FALSE,
            aes(y = outcome,
                x = 0.4,
                hjust = 0,
                family = "RobotoCondensed-Regular",
                label = outcome),
            data = \(x) select(x, outcome) |> unique())


# Substitution ------------------------------------------------------------


diet_scores <- read_csv(here::here("data/MLS_diet_scores_regr_coeff.csv")) |>
  slice(2:5) |>
  rename(beta = estimate, std_error = se, statistic = tval, term = names) |>
  mutate(lcl = beta - (1.96*std_error),
         ucl = beta + (1.96*std_error))





n_nhs <- read_csv(here::here("data/diet_scores_complete_obs_nhs1.csv"))

n_nhs |>
  pivot_longer(everything()) |>
  mutate(name = case_when(name == "N_amed" ~ "aMed",
                          name == "N_ahei" ~ "AHEI",
                          name == "N_alcarb" ~ "ALCD",
                          name == "N_tcarb" ~ "LCD",
                          name == "N_vlcarb" ~ "VLCD"))





  cor_res <- read_csv(here::here("data/diet_score_MLS_cor.csv"))




diet_scores <- cor_res |>
  pivot_longer(-names) |>
  slice(2:6) |>
  mutate(name = case_when(name == "avg_amed" ~ "aMed",
                          name == "avg_ahei" ~ "AHEI",
                          name == "avg_alcarb" ~ "ALCD",
                          name == "avg_tcarb" ~ "LCD",
                          name == "avg_vlcarb" ~ "VLCD")) |>
  ggplot(aes(y = value,
             x = fct_reorder(name, value))) +
  geom_linerange(aes(ymin = 0, ymax = value),
                 lty = 3) +
  geom_hline(yintercept = 0) +
  geom_point(size = 5,
             aes(fill = value),
             shape = 21,
             show.legend = FALSE) +
  scale_fill_gradient2(high = "blue",
                       mid = "white",
                       low = "red",
                       midpoint = 0) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(limits = c(-0.08, 0.15)) +
  labs(y = "Correlation with rMLS",
       x = NULL) +
  theme_light(base_family = "RobotoCondensed-Regular") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(size = 0.1),
        panel.grid = element_line(color = "grey70"),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = "#333333"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_blank())


substitution <-  read_csv(here::here("data/fabian/nutrient_subst_table_pooled.csv")) |>
  slice(2:4)  |>
  mutate(term = case_when(term == "pe_from_ufa_av" ~ " With unsaturated fat",
                          term == "pe_from_prot_av" ~ "With protein",
                          term == "pe_from_carbsav" ~ "With carbohydrates")) |>
  slice(1) |>
  ggplot(aes(y = term,
             x = beta,
             xmin = LCL,
             xmax = UCL)) +
  geom_pointrange(shape = "square",
                  fatten = 8) +
  labs(y = NULL,
       x = "Change in rMLS by substitution of saturated\nwith unsaturated fat (95%-CI)") +
  scale_x_continuous(limits = c(0, 1.15),
                     labels = seq(0, 1.2, 0.2),
                     breaks = seq(0, 1.2, 0.2),
                     expand = c(0, 0)) +
  theme_light(base_family = "RobotoCondensed-Regular") +
  theme(panel.background = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(size = 0.1),
        panel.grid = element_line(color = "grey70"),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = "#333333"),
        panel.grid.major.y = element_blank(),
        axis.ticks = element_blank())


# Multipanel figure -------------------------------------------------------

library(patchwork)


layout <-
  "AACCC
   BBCCC
   BBCCC
   BBCCC"

figure <- substitution + diet_scores + nhs +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = "A")


figure

# Set path
path_to <- ("/home/fabian/alle_shortcut/!MEP/Projekte/EPIC-Potsdam/Diabetes/Lipidomics/MultiLipidScore/AIP/")

path <- str_c(path_to, "Figure_5")
# Save as pdf
ggsave(plot = figure,
       glue::glue("{path}.pdf"),
       device = cairo_pdf,
       width = 25,
       height = 10,
       units = "cm")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)


