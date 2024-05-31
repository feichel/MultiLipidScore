# Functions and packages
library(tidyverse)
library(patchwork)


source(here::here("R/00_functions.R"))

# load weights
comb_weights <- read_rds(here::here("data/weigths_comb_score.rds"))
sep_weights <- read_rds(here::here("data/weigths_separate_scores.rds"))
weigths_all <- bind_rows(comb_weights, sep_weights)

# get clean lipid names
fa_species <- read_csv("/home/fabian/lipidomics.cardiometabolic.diseases/doc/Tables/full_results.csv") %>%
  select(lipid_display = lipid, base, length1, db1, type) %>%
  filter((type == "molecular species" & base %in% c("CE", "Cer", "dhCer", "MG",
                                                    "HexCer", "LacCer", "SM", "LPC",
                                                    "LPE", "FFA") ) |
           type ==  "within-class FA sum",
         !str_detect(lipid_display, "FFA")) %>%
  unique()



# Selecte lipids ---------------------------------------------------------

selected_fig <- fa_species %>%
  mutate(fa = str_c(length1, ":", db1) %>%
           fct_inorder(),
         diet = if_else(lipid_display %in% comb_weights$lipid_display, "Selected", "Not selected")) %>%
  ggplot(aes(y = fct_rev(fa),
             x = reorder_base_new(base))) +
  scale_fill_manual(values = c("Selected" = "#264653",
                               "Not selected" = "grey80")) +
  geom_point(aes(fill = diet),
             shape = 21,
             size = 6,
             show.legend = TRUE) +
  scale_x_discrete(position = "top") +
  labs(x = NULL,
       y = "Class-specific FA coverage",
       fill = NULL) +
  theme_light(base_family = "RobotoCondensed-Regular") +
  theme(panel.grid.major = element_line(color = "#333333",
                                        size = 0.1),
        legend.position = c(0.3, -0.05),
        axis.text.x.top = element_text(angle = 60,
                                       hjust = 0,
                                       vjust = 0),
        axis.text = element_text(color = "#333333"),
        legend.direction = "horizontal",
        legend.key.size = unit(units = "cm", x = 0.2),
        text = element_text(color = "#333333"))




# Weights -----------------------------------------------------------------

order <- weigths_all %>%
  filter(term == "dietUFA") %>%
  select(lipid_display, estimate) %>%
  left_join(fa_species) %>%
  mutate(base = reorder_base_new(base)) %>%
  arrange(fct_rev(base), -estimate)


df_figure_weights <- comb_weights %>%
  select(lipid_display, model, base) %>%
  left_join(fa_species) %>%
  mutate(coeffs = map(model, ~broom::tidy(., conf.int = TRUE))) %>%
  unnest(coeffs) %>%
  filter(term == "dietUFA") %>%
  mutate(weight = sprintf('%.2f', round(estimate, 2)))

write_rds(df_figure_weights, here::here("data/lipids_and_weights.rds"))

estimate_figure <- df_figure_weights %>%
  ggplot(aes(y = fct_relevel(lipid_display, order$lipid_display))) +
  geom_pointrange(aes(x = estimate,
                      xmin = conf.low,
                      xmax = conf.high),
                  shape = "square") +
  labs(y = NULL,
       x = "UFA intervention effect relative\nto SFA-rich diet (95%-CI) [log(μM)]",
       fill = NULL) +
  theme_light(base_family = "RobotoCondensed-Regular") +
  geom_vline(xintercept = 0, lty = 2) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_text(color = "#333333"),
        axis.text.y = element_text(size = 7),
        panel.grid.major.y = element_blank()) +
  ggforestplot::geom_stripes()


weight_column <- df_figure_weights %>%
  ggplot(aes(y = fct_relevel(lipid_display, order$lipid_display))) +
  geom_text(aes(label = weight,
                x = "Weight"),
            size = 3) +
  labs(x = NULL,
       y = NULL) +
  scale_x_discrete(position = "top") +
  theme(axis.text.y = element_blank(),
        axis.text.x.top = element_text(color = "black",
                                       size = 8),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(color = "black"),
        plot.margin = margin(0,0,0,0))




# DIVAS target composition ------------------------------------------------

target_comp <- tribble(
  ~diet, ~sfa_t, ~ufa_t,
  "SFA-rich diet", 17, 15,
  "UFA-rich diet", 9, 23)



target_comp_fig <- target_comp %>%
  pivot_longer(-c(diet)) %>%
  ggplot(aes(y = fct_rev(diet),
             x = value,
             fill = fct_rev(name))) +
  geom_bar(position = "stack", stat = "identity",
           orientation = "y",
           show.legend = FALSE,
           alpha = 0.5) +
  facet_wrap(~diet,
             scales = "free_y",
             ncol = 1) +
  labs(x = "Target % of total energy intake",
       y = NULL,
       fill = NULL) +
  scale_fill_manual(values = c("sfa_t" = "#e9c46a",
                               "ufa_t" = "#2a9d8f")) +
  theme_light(base_family = "RobotoCondensed-Regular") +
  coord_cartesian(expand = FALSE) +
  theme(panel.grid = element_blank(),
        text = element_text(color = "black"),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(color = "black"),
        axis.ticks.y = element_blank()) +
  geom_text(data = tibble(diet = "SFA-rich diet",
                          value = 17/2,
                          name = "sfa_t"),
            label = "SFA",
            color = "black") +
  geom_text(data = tibble(diet = "SFA-rich diet",
                          value = 17 + (15 / 2),
                          name = "sfa_t"),
            label = "UFA",
            color = "black")




# Effects on biomarkers and lipidomics score ------------------------------

combined_score_divas <- read_rds(here::here("data/combined_score_divas.rds")) %>%
  select(volunteer, comb_score_v1, comb_score_v2)


divas_trad <- readxl::read_xlsx(here::here("data/Lipids_DIVAS_for_Fabian.xlsx")) %>%
  janitor::clean_names() %>%
  left_join(combined_score_divas) %>%
  mutate(across(c(total_cholesterol_v1:il6_v2), ~na_if(., 0)),
         across(c(total_cholesterol_v1:il6_v2), ~log(.) %>%
                  scale() %>%
                  as.numeric()),
         across(contains("comb_score"), ~scale(.) %>%
                  as.numeric()),
         diet = if_else(diet == "A", "SFA", "UFA"))


divas_effects_fig <- tibble(var = c("hdl",
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
  mutate(p.value_adj = p.adjust(p.value, method = "fdr")) %>%
  mutate(var = case_when(var == "hdl" ~ "HDL-C",
                         var == "non_hdl" ~ "Non-HDL-C",
                         var == "glucose" ~ "Fasting glucose",
                         var == "crp" ~ "hsCRP",
                         var == "il6" ~ "IL-6",
                         var == "homa_ir" ~ "HOMA-IR",
                         var == "comb_score" ~ "MLS",
                         TRUE ~ "TG")) %>%
  mutate(var = fct_relevel(var, c("HDL-C", "Non-HDL-C", "TG", "Fasting glucose", "HOMA-IR",
                                  "hsCRP", "IL-6", "MLS"))) %>%
  ggplot() +
  ggforestplot::geom_stripes(aes(y = fct_rev(var))) +
  geom_vline(xintercept = 0,
             lty = 2) +
  labs(x = "UFA intervention effect relative\nto SFA-rich diet (95%-CI) [Z-scores]",
       y = NULL,
       fill = NULL) +
  scale_x_continuous(breaks = scales::breaks_pretty()) +
  theme_light(base_family = "RobotoCondensed-Regular") +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  facet_grid(rows = vars(var == "MLS"),
             space = "free",
             scales = "free") +
  geom_pointrange(aes(x = if_else(var == "MLS", estimate * -1, estimate),
                      y = var,
                      xmin = if_else(var == "MLS", conf.low * -1, conf.low),
                      xmax = if_else(var == "MLS", conf.high * -1, conf.high)),
                  orientation = "y",
                  shape = "square",
                  size = 0.5,
                  fatten = 5,
                  position = position_dodge2(width = 0.5),
                  stroke = 0.5) +
  theme(panel.background = element_blank(),
        strip.background.y = element_blank(),
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
        panel.grid.major.y = element_blank(),
        strip.text = element_blank())

layout <- "
AABB
CCBB
CCBB
CCBB
CCDD
CCDD
"

figure <- target_comp_fig +  estimate_figure + selected_fig + divas_effects_fig +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = "A")



path_to <- ("/home/fabian/alle_shortcut/!MEP/Projekte/EPIC-Potsdam/Diabetes/Lipidomics/MultiLipidScore/AIP/")

path <- str_c(path_to, "Figure_2")

# Save as pdf
ggsave(plot = figure,
       glue::glue("{path}.pdf", ),
       device = cairo_pdf,
       width = 22,
       height = 22,
       units = "cm")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)




# extended figure 1
achieved_comp <- tribble(
  ~diet, ~t_a, ~sfa_t, ~ufa_t,
  "SFA-rich diet", "Target", 17, 15,
  "UFA-rich diet", "Target", 9, 23,
  "SFA-rich diet", "Achieved", 17.6, 14.5,
  "UFA-rich diet", "Achieved", 8.05, 24)


figure <- achieved_comp |>
  pivot_longer(-c(diet, t_a)) %>%
  ggplot(aes(y = t_a,
             x = value,
             fill = fct_rev(name))) +
  geom_bar(position = "stack", stat = "identity",
           orientation = "y",
           show.legend = FALSE,
           alpha = 0.5) +
  facet_wrap(~diet,
             scales = "free_y",
             ncol = 1) +
  labs(x = "Target % of total energy intake",
       y = NULL,
       fill = NULL) +
  scale_fill_manual(values = c("sfa_t" = "#e9c46a",
                               "ufa_t" = "#2a9d8f")) +
  theme_light(base_family = "RobotoCondensed-Regular") +
  coord_cartesian(expand = FALSE) +
  theme(panel.grid = element_blank(),
        text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(color = "black"),
        axis.ticks.y = element_blank()) +
  geom_text(data = tibble(diet = "SFA-rich diet",
                          value = 17/2,
                          t_a = "Target",
                          name = "sfa_t"),
            label = "SFA",
            color = "black") +
  geom_text(data = tibble(diet = "SFA-rich diet",
                          value = 17 + (15 / 2),
                          name = "sfa_t",
                          t_a = "Target"),
            label = "UFA",
            color = "black")



path_to <- ("/home/fabian/alle_shortcut/!MEP/Projekte/EPIC-Potsdam/Diabetes/Lipidomics/MultiLipidScore/AIP/")

path <- str_c(path_to, "Extended_Figure_1")

# Save as pdf
ggsave(plot = figure,
       glue::glue("{path}.pdf", ),
       device = cairo_pdf,
       width = 15,
       height = 10,
       units = "cm")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)
