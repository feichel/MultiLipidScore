library(tidyverse)


data_for_analysis <- read_rds(here::here("data/data_for_analysis.rds"))

combined_score <- read_rds(here::here("data/predimed_score.rds")) %>%
  rename(predimed_score = comb_score)


panel_a <- data_for_analysis %>%
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
        axis.line = element_line()) +
  labs(y = "rMLS (Z-scores)",
       x = "MLS (Z-scores)") +
  ggpubr::stat_cor(p.accuracy = 0.001,
                   label.x = 1.6,
                   label.y = -5,
                   method = "spearman",
                   cor.coef.name = "rho")

temp <- data_for_analysis %>%
  select(omics_id, comb_score) %>%
  left_join(combined_score)


cor.test(temp$comb_score, temp$predimed_score,
         method = "spearman")



df <- data_for_analysis %>%
  select(omics_id, comb_score) %>%
  left_join(combined_score) |>
  mutate(
    across(
      c(comb_score, predimed_score), ~ scale(.) |>
        as.numeric()
    ),
    means = (comb_score + predimed_score) / 2,
    diffs = comb_score - predimed_score
  )

sd_diffs <- sd(df$diffs)
bias <- mean(df$diffs)



# Calculating limits of agreement
upper_loa <- bias + 2 * sd_diffs
lower_loa <- bias - 2 * sd_diffs

# Plotting using ggplot2
panel_b <- ggplot(df, aes(x = means, y = diffs)) +
  geom_point() +
  labs(y = "Difference",
       x = "Mean") +
  theme_light(base_family = "RobotoCondensed-Regular") +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()) +
  geom_hline(yintercept = bias, color = "red") +
  geom_hline(yintercept = upper_loa, color = "black", linetype = "dashed") +
  geom_hline(yintercept = lower_loa, color = "black", linetype = "dashed") +
  expand_limits(x = c(min(df$means), max(df$means) * 1.1)) +
  annotate("text", x = Inf, y = upper_loa + 0.1, label = "+2×SD", hjust = 1.1) +
  annotate("text", x = Inf, y = lower_loa + 0.1, label = "-2×SD", hjust = 1.1)

library(patchwork)

figure <- panel_a + panel_b + plot_annotation(tag_levels = "A")


  # Set path
  path_to <- ("/home/fabian/alle_shortcut/!MEP/Projekte/EPIC-Potsdam/Diabetes/Lipidomics/MultiLipidScore/AIP/")

path <- str_c(path_to, "Extended_Figure_5")

# Save as pdf
ggsave(plot = figure,
       glue::glue("{path}.pdf"),
       device = cairo_pdf,
       width = 30,
       height = 15,
       units = "cm")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)
