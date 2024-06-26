library(tidyverse)



# Predimed figure ---------------------------------------------------------
# input tables provided

stratum <- 0.06
n_total <- 0.2
n_cases <- 0.33
hr_text <- 2.6
y <- 3.7

all <- read_csv(here::here("data/Intervention_effects_according_to_lip_score_1.csv"),
                     name_repair = janitor::make_clean_names) %>%
  mutate(strata = case_when(strata == "BEN_MLS" ~ "Beneficial",
                            strata == "ADV_MLS" ~ "Adverse",
                            TRUE ~ "Pooled") %>%
           fct_relevel(c("Pooled","Adverse", "Beneficial"))) %>%
  ggplot() +
  theme_light(base_family = "RobotoCondensed-Regular") +
  labs(y = NULL,
       x = NULL) +
  ggforestplot::geom_stripes(aes(y = strata),
                             show.legend = FALSE) +
  scale_y_discrete(expand = c(0,0)) +
  geom_pointrange(aes(y = strata,
                      x = hr,
                      xmin = lcl,
                      xmax = ucl),
                  shape = "square") +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank()) +
  coord_cartesian(clip = "off") +
  scale_x_log10(breaks = 1, limits = c(0.06, 4)) +
  geom_text(inherit.aes = FALSE,
            aes(y = strata,
                x = stratum,
                hjust = 0,
                family = "RobotoCondensed-Regular",
                label = strata)) +
  geom_text(inherit.aes = FALSE,
            aes(y = strata,
                x = n_total,
                family = "RobotoCondensed-Regular",
                label = n)) +
  geom_text(inherit.aes = FALSE,
            aes(y = strata,
                x = n_cases,
                family = "RobotoCondensed-Regular",
                label = nevent)) +
  geom_text(inherit.aes = FALSE,
            aes(y = strata,
                x = hr_text,
                family = "RobotoCondensed-Regular",
                label = hr_95_percent_ci)) +
  annotate(geom = "text",
           y = y,
           fontface = 2,
           hjust = 0,
           x = stratum,
           family = "RobotoCondensed-Regular",
           label = "rMLS stratum") +
  annotate(geom = "text",
           y = y,
           x = n_total,
           family = "RobotoCondensed-Regular",
           fontface = 2,
           label = "Total n") +
  annotate(geom = "text",
           y = y,
           x = n_cases,
           family = "RobotoCondensed-Regular",
           fontface = 2,
           label = "Cases") +
  annotate(geom = "text",
           y = y,
           x = hr_text,
           family = "RobotoCondensed-Regular",
           fontface = 2,
           label = "HR (95%-CI)") +
  annotate(geom = "segment",
           y = 0.2,
           yend = 0.2,
           x = 0.9,
           xend = 0.7,
           arrow = arrow(type = "closed",
                         length = unit(0.03, "npc"))) +
  annotate(geom = "segment",
           y = 0.2,
           yend = 0.2,
           x = 1.1,
           xend = 1.4,
           arrow = arrow(type = "closed",
                         length = unit(0.03, "npc"))) +
  annotate(geom = "text",
           y = -0.1,
           x = 0.9,
           hjust = 1,
           family = "RobotoCondensed-Regular",
           label = "Favors intervention") +
  annotate(geom = "text",
           y = -0.1,
           x = 1.1,
           hjust = 0,
           family = "RobotoCondensed-Regular",
           label = "Favors control")





stratum_2 <- 0.06
hr_text_2 <- 0.15
y_2 <- 3.7

by_arm <- read_csv(here::here("data/MD_T2D_effects_according_to_lip_score.csv"),
                   name_repair = janitor::make_clean_names) |>
  mutate(strata = case_when(strata == "HIGH_LIPS" ~ "Beneficial",
                            strata == "LOW_LIPS" ~ "Adverse",
                            TRUE ~ "Pooled") |>
           fct_relevel(c("Pooled","Adverse", "Beneficial")),
         group = if_else(str_detect(group, "EVOO"), "EVOO", "Nuts")) |>
  ggplot() +
  theme_light(base_family = "RobotoCondensed-Regular") +
  labs(y = NULL,
       x = NULL,
       fill = "Mediterranean diet intervention") +
  ggforestplot::geom_stripes(aes(y = strata),
                             show.legend = FALSE) +
  ggsci::scale_fill_d3() +
  scale_y_discrete(expand = c(0,0)) +
  geom_pointrange(aes(y = strata,
                      x = hr,
                      xmin = lcl,
                      xmax = ucl,
                      fill = group),
                  shape = 21,
                  position = position_dodge(width = 0.5)) +
  theme(panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        legend.position = "bottom") +
  coord_cartesian(clip = "off") +
  scale_x_log10(breaks = 1, limits = c(0.06, 2)) +
  geom_text(inherit.aes = FALSE,
            aes(y = strata,
                x = stratum_2,
                hjust = 0,
                family = "RobotoCondensed-Regular",
                label = strata)) +
  geom_text(inherit.aes = FALSE,
            aes(y = strata,
                x = hr_text_2,
                hjust = 0,
                family = "RobotoCondensed-Regular",
                label = hr_95_percent_ci,
                group = group),
            position = position_dodge(width = 0.5)) +
  annotate(geom = "text",
           y = y_2,
           fontface = 2,
           hjust = 0,
           x = stratum,
           family = "RobotoCondensed-Regular",
           label = "rMLS stratum") +
  annotate(geom = "segment",
           y = 0.2,
           yend = 0.2,
           x = 0.9,
           xend = 0.7,
           arrow = arrow(type = "closed",
                         length = unit(0.03, "npc"))) +
  annotate(geom = "segment",
           y = 0.2,
           yend = 0.2,
           x = 1.1,
           xend = 1.4,
           arrow = arrow(type = "closed",
                         length = unit(0.03, "npc"))) +
  annotate(geom = "text",
           y = -0.1,
           x = 0.9,
           hjust = 1,
           family = "RobotoCondensed-Regular",
           label = "Favors intervention") +
  annotate(geom = "text",
           y = -0.1,
           x = 1.1,
           hjust = 0,
           family = "RobotoCondensed-Regular",
           label = "Favors control")





library(patchwork)

figure <- all / by_arm + patchwork::plot_annotation(tag_levels = "A")



# Set path
path_to <- ("/home/fabian/alle_shortcut/!MEP/Projekte/EPIC-Potsdam/Diabetes/Lipidomics/MultiLipidScore/AIP/")

path <- str_c(path_to, "Figure_6")
# Save as pdf
ggsave(plot = figure,
       glue::glue("{path}.pdf"),
       device = cairo_pdf,
       width = 25,
       height = 16,
       units = "cm")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)


