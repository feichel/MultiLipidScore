library(tidyverse)


data_for_analysis <- read_rds(here::here("data/data_for_analysis.rds"))

combined_score <- read_rds(here::here("data/predimed_score.rds")) %>%
  rename(predimed_score = comb_score)


figure <- data_for_analysis %>%
  select(omics_id, comb_score) %>%
  left_join(combined_score) %>%
  ggplot(aes(x = scale(comb_score),
             y = scale(predimed_score))) +
  geom_point() +
  geom_abline(color = "red") +
  coord_equal(xlim = c(-5.5, 5.5),
              ylim = c(-5.5, 5.5),
              expand = FALSE) +
  scale_y_continuous(breaks = seq(-5, 5, 2.5)) +
  scale_x_continuous(breaks = seq(-5, 5, 2.5)) +
  theme_light(base_family = "RobotoCondensed-Regular") +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(),
        plot.margin = unit(c(1,1,1,1), "lines")) +
  labs(y = "rMLS (Z-scores)",
       x = "MLS (Z-scores)") +
  ggpubr::stat_cor(p.accuracy = 0.001,
                   label.x = 1.6,
                   label.y = -5)


# Set path
path <- here::here("doc", "Figures", "ext_fig1")

# Save as pdf
ggsave(plot = figure,
       glue::glue("{path}.pdf"),
       device = cairo_pdf,
       width = 15,
       height = 15,
       units = "cm")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)
