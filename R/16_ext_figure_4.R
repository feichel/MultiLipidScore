library(tidyverse)

# Substitution extended ---------------------------------------------------

paths <- c(here::here("data/fabian/nutrient_subst_table_pooled.csv"),
           here::here("data/fabian/nutrient_subst_table_nhs1.csv"),
           here::here("data/fabian/nutrient_subst_table_nhs2.csv")) |>
  set_names(c("Pooled", "NHS1", "NHS2"))

extended_substitution <- map_dfr(paths, read_csv,
                                 .id = "path") |>
  filter(str_detect(term, "^pe_from")) |>
  mutate(term = case_when(term == "pe_from_ufa_av" ~ " For unsaturated fat",
                          term == "pe_from_prot_av" ~ "For protein",
                          term == "pe_from_carbsav" ~ "For carbohydrates")) %>%
  ggplot(aes(y = term,
             x = beta,
             xmin = LCL,
             xmax = UCL,
             fill = fct_rev(path))) +
  geom_pointrange(shape = 22,
                  fatten = 8,
                  position = position_dodge(width = 0.5)) +
  labs(y = NULL,
       fill = NULL,
       x = "Change in rMLS by substitution of saturated fat (95%-CI)") +
  scale_x_continuous(limits = c(0, 1.15),
                     labels = seq(0, 1.2, 0.2),
                     breaks = seq(0, 1.2, 0.2),
                     expand = c(0, 0)) +
  scale_fill_manual(breaks = c("NHS1", "NHS2", "Pooled"),
                    values = c("yellow", "red", "blue3")) +
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


# Set path
path <- here::here("doc", "Figures", "predimed_extended")

# Save as pdf
ggsave(plot = extended_substitution,
       glue::glue("{path}.pdf"),
       device = cairo_pdf,
       width = 15,
       height = 10,
       units = "cm")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)



