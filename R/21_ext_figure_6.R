library(tidyverse)


nhs <- read_csv(here::here("data/fa_intake_dist_nhs.csv"))
nhs2 <- read_csv(here::here("data/fa_intake_dist_nhs2.csv"))
read_csv(here::here("data/fa_intake_dist_nhs_summary_table.csv"))
read_csv(here::here("data/fa_intake_dist_nhs2_summary_table.csv"))

bind_rows(nhs, nhs2) |>
  na.omit() |>
  nrow()
# [1] 10381


temp_figure <-bind_rows(nhs, nhs2) |>
  na.omit() |>
  ggplot(aes(x = FATpe,
             y = UFAtoSFAratio)) +
  geom_point(fill = "#2c73d2",
             shape = 21) +
  geom_hline(yintercept = 0.8823529, lty = 2) +
  geom_hline(yintercept = 2.555556, lty = 2) +
  geom_vline(xintercept = 36, lty = 2) +
  annotate("segment",
           x = 36,
           xend = 36,
           y = 0.8823529,
           yend = 2.555556 ,
           alpha = .9,
           size = 2,
           color = "red",
           arrow = arrow(ends = "both", angle = 90)) +
  annotate("label",
           label = "DIVAS target difference\nin UFA-to-SFA ratio",
           size = 5,
           x = 42,
           y = (0.8823529 + 2.555556) /2,
           color = "red") +
  annotate("label",
           label = "SFA-rich diet target",
           size = 5,
           x = 10,
           hjust = "left",
           y = (0.8823529),
           color = "black") +
  annotate("label",
           label = "UFA-rich diet target",
           size = 5,
           x = 10,
           hjust = "left",
           y = (2.555556),
           color = "black") +
  annotate("label",
           label = "Total fat target",
           size = 5,
           x = 36,
           hjust = "left",
           y = 0.4,
           color = "black") +
  scale_x_continuous(expand = c(0, 0.5)) +
  labs(y = "Dietary UFA-to-SFA ratio",
       x = "Total energy from fat [E%]") +
  theme_classic(base_family = "RobotoCondensed-Regular") +
  theme(panel.background = element_blank(),
        legend.position = "bottom",
        axis.ticks.y = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = "#333333"),
        strip.background = element_rect(fill = "#333333"),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_text(vjust = 0))


figure <- ggExtra::ggMarginal(temp_figure,
                    type = "histogram",
                    size = 5 ,
                    fill = "#ff8066",
                    xparams = list(bins = 100),
                    yparams = list(bins = 100))

# Set path
path_to <- ("/home/fabian/alle_shortcut/!MEP/Projekte/EPIC-Potsdam/Diabetes/Lipidomics/MultiLipidScore/AIP/")

path <- str_c(path_to, "Extended_Figure_6")

# Save as pdf
ggsave(plot = figure,
       glue::glue("{path}.pdf"),
       device = cairo_pdf,
       width = 20,
       height = 20,
       units = "cm")


# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)
