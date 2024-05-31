


source(here::here("R/00_functions.R"))
library(tidyverse)




display_names <- read_rds(here::here("data/weigths_comb_score.rds")) %>%
  select(lipid, lipid_display)


single_lip_res <- read_rds(here::here("data/cluster.info.rds"))



order <- single_lip_res %>%
  select(lipid, lipid_display, louvain) %>%
  separate(lipid_display, into = c("base", "l", "db"),
           remove = FALSE,
           extra = "drop") %>%
  mutate(l = str_remove(l, "FA"),
         base = reorder_base_new(base)) %>%
  arrange(louvain, fct_rev(base), l, db)



netcoupler_direct_effects <- read_rds(file = here::here("doc/netcoupler_direct_effects.rds"))


single_assocs <- single_lip_res %>%
  mutate(nc = case_when(outcome == "CVD" & lipid %in% netcoupler_direct_effects$cvd ~ TRUE,
                        outcome == "Type 2 diabetes" & lipid %in% netcoupler_direct_effects$diab ~ TRUE,
                        TRUE ~ FALSE)) %>%
  ggplot() +
  facet_wrap(~outcome,
             ncol = 2) +
  geom_vline(xintercept = 1,
             lty = 2) +
  labs(x = "HR per SD (95%-CI)",
       y = NULL,
       fill = "Cluster") +
  theme_light(base_family = "RobotoCondensed-Regular") +
  scale_fill_manual(values = paletteer::paletteer_d("fishualize::Cirrhilabrus_solorensis",
  direction = 1)) +
  theme(panel.background = element_blank(),
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
        strip.background.y = element_blank(),
        strip.text.y = element_blank()) +
  geom_pointrange(aes(x = estimate,
                      y = fct_rev(fct_relevel(lipid_display, unique(order$lipid_display))),
                      xmin = conf.low,
                      xmax = conf.high,
                      fill = louvain),
                  shape = 22,
                  show.legend = TRUE) +
  geom_point(data = function(x) filter(x, nc == TRUE),
             shape = 8,
             color = "red",
             aes(x = conf.high + (conf.high * 0.2),
                 y = fct_rev(fct_relevel(lipid_display, unique(order$lipid_display))))) +
  scale_x_log10(breaks = c(0.3, 1, 3, 8))



# set path
path_to <- ("/home/fabian/alle_shortcut/!MEP/Projekte/EPIC-Potsdam/Diabetes/Lipidomics/MultiLipidScore/AIP/")

path <- str_c(path_to, "Extended_Figure_10")


# Save as pdf
ggsave(plot = single_assocs,
       glue::glue("{path}.pdf"),
       device = cairo_pdf,
       width = 15,
       height = 20,
       units = "cm")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)
