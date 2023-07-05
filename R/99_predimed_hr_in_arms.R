library(tidyverse)


figure <- read_csv(here::here("data/LS_T2D_by_MDintervention.csv"),
                name_repair = janitor::make_clean_names) %>%
  ggplot() +
  ggforestplot::geom_stripes(aes(y = fct_rev(strata)),
                             show.legend = FALSE) +
  scale_y_discrete(expand = c(0, 0.5)) +
  geom_vline(xintercept = 1,
             lty = 2) +
  geom_pointrange(aes(x = hr,
                      y = fct_rev(strata),
                      xmin = lcl,
                      xmax = ucl),
                  orientation = "y",
                  shape = "square") +
  labs(x = "Hazard ratio per SD (95% CI)",
       y = NULL) +
  theme_light() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1),
        panel.grid = element_line(color = "grey70"),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = "#333333"),
        strip.background = element_rect(fill = "#333333"),
        panel.grid.major.y = element_blank()) +
  scale_x_log10(breaks = c(seq(0.6, 1, 0.1),
                           map_dbl(seq(0.5, 1, 0.1), ~1/.) %>% round(digits = 1)),
                expand = c(0, 0.01)) +
  geom_text(inherit.aes = FALSE,
            aes(y = fct_rev(strata),
                x = 0.4,
                label = hr_95_percent_ci),
            hjust = 0)


# Set path
path <- here::here("doc", "Figures", "predimed_supplement")

# Save as pdf
ggsave(plot = figure,
       glue::glue("{path}.pdf"),
       device = cairo_pdf,
       width = 22,
       height = 7,
       units = "cm")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)

