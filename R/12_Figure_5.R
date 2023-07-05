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
  baseline_stroke_adj
) %>%
  janitor::clean_names() %>%
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
        legend.position = "bottom") +
  labs(y = NULL,
       fill = NULL,
       x = "Relative risk per SD increase in rMLS (95%-CI)") +
  ggsci::scale_fill_npg(breaks = c("BMI-unadjusted", "BMI-adjusted")) +
  scale_x_log10(breaks = c(0.6, 0.7,0.8, 0.9, 1), limits = c(0.4, 1.05)) +
  geom_text(inherit.aes = FALSE,
            aes(y = outcome,
                x = 0.4,
                hjust = 0,
                family = "RobotoCondensed-Regular",
                label = outcome),
            data = \(x) select(x, outcome) |> unique())

# not used
# bind_rows(
#   change_t2d,
#   baseline_t2d,
#   baseline_stroke,
#   baseline_stroke_adj
# ) %>%
#   janitor::clean_names() %>%
#   separate(term, into = c("type", "outcome"), remove = FALSE) %>%
#   mutate(type = str_to_sentence(type),
#          outcome = str_replace(outcome, "t2d", "Type 2 diabetes") %>%
#            str_replace("stroke", "Stroke")) %>%
#   ggplot() +
#   geom_vline(xintercept = 1) +
#   geom_pointrange(aes(y = term,
#                       x = hr,
#                       xmin = lcl,
#                       xmax = ucl),
#                   shape = "square") +
#   ggforestplot::geom_stripes(aes(y = term),
#                              show.legend = FALSE) +
#   theme_light(base_family = "RobotoCondensed-Regular") +
#   theme(panel.background = element_blank(),
#         panel.border = element_blank(),
#         axis.text.y = element_blank(),
#         panel.grid.major = element_line(size = 0.1),
#         panel.grid = element_line(color = "grey70"),
#         panel.grid.minor = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.x = element_text(color = "#333333"),
#         panel.grid.major.y = element_blank()) +
#   labs(y = NULL,
#        x = "Hazard ratio per SD increase (95%-CI)") +
#   scale_x_log10(breaks = c(0.6, 0.7,0.8, 0.9, 1), limits = c(0.4, 1)) +
#   geom_text(inherit.aes = FALSE,
#             aes(y = term,
#                 x = 0.4,
#                 hjust = 0,
#                 family = "RobotoCondensed-Regular",
#                 label = outcome)) +
#   geom_text(inherit.aes = FALSE,
#             aes(y = term,
#                 x = 0.5,
#                 hjust = 0,
#                 family = "RobotoCondensed-Regular",
#                 label = type)) +
#   annotate(geom = "text",
#            y = y,
#            fontface = 2,
#            hjust = 0,
#            x = 0.5,
#            family = "RobotoCondensed-Regular",
#            label = "Variable") +
#   annotate(geom = "text",
#            y = y,
#            x = 0.4,
#            hjust = 0,
#            family = "RobotoCondensed-Regular",
#            fontface = 2,
#            label = "Outcome") +
#   coord_cartesian(clip = "off")


# Substitution ------------------------------------------------------------


substitution <-  read_csv(here::here("data/fabian/nutrient_subst_table_pooled.csv")) %>%
  slice(2:4)  %>%
  mutate(term = case_when(term == "pe_from_ufa_av" ~ " For unsaturated fat",
                          term == "pe_from_prot_av" ~ "For protein",
                          term == "pe_from_carbsav" ~ "For carbohydrates")) %>%
  ggplot(aes(y = term,
             x = beta,
             xmin = LCL,
             xmax = UCL)) +
  geom_pointrange(shape = "square",
                  fatten = 8) +
  labs(y = NULL,
       x = "Change in rMLS by substitution\nof saturated fat (95%-CI)") +
  scale_x_continuous(limits = c(0, 1.15),
                     labels = seq(0, 1.2, 0.2),
                     breaks = seq(0, 1.2, 0.2),
                     expand = c(0, 0)) +
  theme_light(base_family = "RobotoCondensed-Regular") +
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(size = 0.1),
        panel.grid = element_line(color = "grey70"),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = "#333333"),
        panel.grid.major.y = element_blank())




# Multipanel figure -------------------------------------------------------

library(patchwork)


layout <- "
AABBB
"

figure <- substitution +
  nhs +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = "A")



# Set path
path <- here::here("doc", "Figures", "Figure_5")

# Save as pdf
ggsave(plot = figure,
       glue::glue("{path}.pdf"),
       device = cairo_pdf,
       width = 22,
       height = 10,
       units = "cm")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)


