library(tidyverse)



## data
sls <- read_rds(here::here("data/sphingolipid_score_epic.rds")) %>%
  filter(ia_subcohort == 1) |>
  select(omics_id, sphingolipid_score)


mls <- read_rds(here::here("data/data_for_analysis.rds")) %>%
  filter(ia_subcohort == 1) |>
  select(omics_id, comb_score)


cor.test(sls$sphingolipid_score, mls$comb_score,
         method = "spearman")


panel_b <- sls %>%
  left_join(mls) %>%
  ggplot(aes(x = scale(comb_score),
             y = scale(sphingolipid_score))) +
  geom_point() +
  geom_abline(color = "red") +
  scale_y_continuous(breaks = seq(-5, 5, 2.5), limits = c(-5.5, 5.5)) +
  scale_x_continuous(breaks = seq(-5, 5, 2.5), limits = c(-5.5, 5.5)) +
  theme_light(base_family = "RobotoCondensed-Regular") +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()) +
  labs(y = "Sphingolipid-Score (Z-score)",
       x = "MLS (Z-score)") +
  ggpubr::stat_cor(p.accuracy = 0.001,
                   label.x = 0.5,
                   label.y = -5,
                   method = "spearman",
                   cor.coef.name = "rho")



lipid_results <- read_rds(here::here("data/lipid_results.rds")) |>
  mutate(lipid = recode(lipid,
                        `cer_fa18_0` = "Cer(18:0)",
                        `dicer_fa18_0` = "dhCer(18:0)",
                        `dicer_fa20_0` = "dhCer(20:0)",
                        `dicer_fa22_0` = "dhCer(22:0)",
                        `dicer_fa24_1` = "dhCer(24:1)",
                        `glucer_fa18_0` = "HexCer(18:0)",
                        `sm_fa14_0` = "SM(14:0)"))

comb_weights <- read_rds(here::here("data/weigths_z_comb_score.rds")) |>
  select(lipid = lipid_display, estimate, std.error) |>
  mutate(study = "DIVAS") |>
  filter(lipid %in% lipid_results$lipid)

lipid_compare <- lipid_results |>
  unnest(coeffs_z) |>
  filter(term == "diet1",
         lipid != "apob") |>
  select(lipid, estimate, std.error) |>
  mutate(study = "Lipogain-2") |>
  bind_rows(comb_weights)


panel_a <- lipid_compare |>
  ggplot(aes(y = lipid)) +
  geom_pointrange(aes(x = estimate,
                      xmin = estimate - 1.96 * std.error,
                      xmax = estimate + 1.96 * std.error,
                      fill = study),
                  shape = 22,
                  fatten = 8,
                  position = position_dodge(width = 0.5)) +
  labs(y = NULL,
       x = "Intervention effect relative\nto control diet (95%-CI) [Z-score]",
       fill = NULL) +
  theme_light(base_family = "RobotoCondensed-Regular") +
  geom_vline(xintercept = 0, lty = 2) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_text(color = "#333333"),
        axis.text.y = element_text(size = 7),
        panel.grid.major.y = element_blank()) +
  ggforestplot::geom_stripes() +
  ggsci::scale_fill_npg()


score_results <- read_rds(here::here("data/score_results.rds")) |>
  filter(model %in% c("sl_1", "sl_2"))


panel_c <- score_results |>
  unnest(coeffs) |>
  filter(term == "diet1",
         !str_detect(model, "_z")) |>
  mutate(adj = if_else(model == "sl_1", "unadjusted", "ApoB-adjusted"),
         lipid = "SL-Score") |>
  select(lipid, adj, estimate, std.error) |>
  ggplot(aes(y = lipid)) +
  geom_pointrange(aes(x = estimate,
                      xmin = estimate - 1.96 * std.error,
                      xmax = estimate + 1.96 * std.error,
                      fill = adj),
                  shape = 22,
                  fatten = 8,
                  position = position_dodge(width = 0.5)) +
  labs(y = NULL,
       x = "Intervention effect on Sphingolipid-Score\nrelative to control (95%-CI)",
       fill = NULL) +
  theme_light(base_family = "RobotoCondensed-Regular") +
  geom_vline(xintercept = 0, lty = 2) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_text(color = "#333333"),
        axis.text.y =element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank()) +
  ggsci::scale_fill_d3()




score_residuals <- score_results |>
  filter(model == "sl_1") |>
  unnest(residuals) |>
  add_column(id = rep(1:60)) |>
  select(id, value = residuals) |>
  mutate(id = if_else(id < 27, id, id + 1),
         name = "sl_score") |>
  add_row(id = 27)


lipid_residuals <- lipid_results |>
  select(lipid, residuals_z) |>
  filter(lipid != "apob") |>
  unnest(residuals_z) |>
  add_column(id = rep(1:60, 7)) |>
  pivot_wider(names_from = lipid,
              values_from = residuals_z) |>
  mutate(id = if_else(id < 27, id, id + 1)) |>
  add_row(id = 27) |>
  pivot_longer(-id)

apob_residuals <- lipid_results |>
  select(lipid, residuals_z) |>
  filter(lipid == "apob") |>
  unnest(residuals_z) |>
  add_column(id = 1:61)





panel_d <- bind_rows(score_residuals) |>
  left_join(apob_residuals) |>
  filter(!is.na(name)) |>
  ggplot(aes(y = residuals_z*-1,
             x = value)) +
  geom_point() +
  geom_abline(color = "red") +
  # coord_equal() +
  scale_x_continuous(breaks = -2:2, limits = c(-2.1, 2.1)) +
  scale_y_continuous(breaks = -2:2, limits = c(-2.1, 2.1)) +
  labs(x = "Sphingolipid-Score (Z-Score)",
       y = "Reduction in ApoB concentration\n(log-transformed, Z-Score)") +
  theme_light(base_family = "RobotoCondensed-Regular") +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line()) +
  ggpubr::stat_cor(p.accuracy = 0.001,
                   label.x = 0.1,
                   label.y = -2,
                   method = "spearman",
                   cor.coef.name = "rho")


temp <- bind_rows(score_residuals) |>
  left_join(apob_residuals) |>
  filter(!is.na(name))

cor.test(temp$residuals_z, temp$value,
         method = "spearman")


library(patchwork)

layout <- c("AB
            AB
            AB
            AB
            AB
            AD
            AD
            CD")




figure <-  panel_a + panel_b + panel_c + panel_d + plot_layout(design = layout) +
  plot_annotation(tag_levels = "A")

figure
# Set path
path_to <- ("/home/fabian/alle_shortcut/!MEP/Projekte/EPIC-Potsdam/Diabetes/Lipidomics/MultiLipidScore/AIP/")

path <- str_c(path_to, "Extended_Figure_3")

# Save as pdf
ggsave(plot = figure,
       glue::glue("{path}.pdf"),
       device = cairo_pdf,
       height = 25,
       width = 25,
       units = "cm",
       bg = "white")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)




