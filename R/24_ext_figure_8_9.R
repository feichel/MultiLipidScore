
# Install netcoupler
# pak::pkg_install("netCoupler/netCoupler")


# Load libraries
library(tidyverse)
library(NetCoupler)
library(ggraph)
library(tidygraph)
library(survival)

# Source functions
source(here::here("R/00_functions.R"))

# Estimate network and identify clusters of lipids
set.seed(123)
network <- read_rds(here::here("data/data_for_clusters.rds")) %>%
  filter(ia_subcohort == 1) %>%
  select(-c(omics_id, ia_subcohort)) %>%
  NetCoupler::nc_estimate_network(everything()) %>%
  mutate(louvain = as.factor(tidygraph::group_louvain()))

write_rds(network, file = here::here("data/network.rds"))



# rename lipids and wrangle network data
lipid_names <- network %>%
  pull(1) %>%
  tibble(lipid = .) |>
  epic_lip_names()



# use SNAhelper to manually move nodes to make it look pretty
# save node positions and hard code

x <- c(1.4368, 0.7345, 1.3802, -0.2686, 0.3683, -0.8742, -0.1354, -1.0144, 0.567, 1.2219, 0.5282, 2.1495, 2.8793,
       2.1654, 2.7751, 2.5173, -0.2004, -0.8838, -0.0808, 0.1981, -1.2998, -0.1465, -0.7681, 1.1789, 0.339, 0.0098,
       0.8929, -0.8732, -1.5003, 2.0517, 0.567, -0.5635, 1.366, -1.6827, -0.7941, -1.6617, -2.0137, -2.2416, -2.0196,
       -2.4879, -2.3864, -1.6512, -0.3561, 1.8368, 1.3699)
y <- c(0.3328, 1.6106, 1.67, 2.2269, 0.6751, 2.4622, -0.3412, -1.0879, -1.294, 3.0399, 0.1869, -2.2125, -0.5309,
       -0.5901, -1.7455, 0.5571, -1.1033, -2.046, -1.7826, -1.1416, -0.2426, 0.708, 0.0691, -0.8513, -2.2166, 2.7817,
       0.8507, 0.6668, -1.2898, 1.8177, 2.0596, -1.3512, -0.4136, -0.4366, 1.6185, 1.4814, -0.2714, 0.1616, 0.9824,
       -0.1558, 0.8401, 0.8047, 1.3241, -0.98, -2.0481)


# plot network with cluster identification
network_figure <- network %>%
  left_join(lipid_names, by = c("name" = "lipid")) %>%
  ggraph(layout = "manual", x = x, y = y) +
  ggraph::geom_edge_link2(width = 0.5,
                          aes(alpha = abs(weight)),
                          end_cap = circle(4.5, 'mm'),
                          start_cap = circle(4.5, 'mm'),
                          show.legend = FALSE) +
  # ggraph::geom_edge_link2(aes(filter = abs(weight) > 0.3,
  #                             label = sprintf('%.2f', round(weight, 2)),
  #                             alpha = abs(weight)),
  #                         label_dodge = unit(2, 'mm'),
  #                         angle_calc = "along",
  #                         label_size = 3,
  #                         end_cap = circle(4.5, 'mm'),
  #                         start_cap = circle(4.5, 'mm')) +
  geom_node_point(aes(colour = louvain),
                  fill = "#333333",
                  shape = 21,
                  size = 10,
                  stroke = 2) +
  theme_graph() +
  theme(legend.position = "bottom") +
  labs(color = "Cluster") +
  paletteer::scale_color_paletteer_d("fishualize::Cirrhilabrus_solorensis") +
  scale_edge_alpha_continuous() +
  geom_node_text(aes(label = base),
                 size = 2,
                 color = "white",
                 position = position_nudge(y = 0.06),
                 fontface = "bold") +
  geom_node_text(aes(label = fa),
                 size = 2,
                 color = "white",
                 position = position_nudge(y = -0.06),
                 fontface = "bold") +
  guides(color = guide_legend(override.aes = list(size = 5)))


# Set path
path_to <- ("/home/fabian/alle_shortcut/!MEP/Projekte/EPIC-Potsdam/Diabetes/Lipidomics/MultiLipidScore/AIP/")

path <- str_c(path_to, "Extended_Figure_8")

# Save as pdf
ggsave(plot = network_figure,
       glue::glue("{path}.pdf"),
       device = cairo_pdf,
       width = 20,
       height = 12,
       units = "cm")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)

# generate cluster-specific scores ---------------------------------------------------

df_lips <- read_rds(here::here("data/data_for_clusters.rds"))

cluster_vector <- paste0("cluster_", seq(1, 5))


#' Function to calculate cluster-specific scores i.e. generate score only from lipids
#' that are in one cluster
#'
#' @param .cluster
#'
#' @return
#' @export
#'
#' @examples
make_cluster_score <- function(.cluster) {

  cluster_lipids <- network %>%
    tidygraph::as_tibble() %>%
    filter(louvain == .cluster) %>%
    pull(name)

  weights <- read_rds(here::here("data/weigths_comb_score.rds")) %>%
    filter(lipid %in% cluster_lipids)

  lipids_to_score <- df_lips %>%
    select(weights$lipid)

  # Calculate score by multiplying weights with respective columns
  return(as.numeric(as.matrix(lipids_to_score) %*% weights$estimate))

}

new_scores <- seq(1, 5) %>%
  set_names(cluster_vector) %>%
  map_dfc(make_cluster_score) %>%
  bind_cols(omics_id = df_lips$omics_id, .)



pheno_df <- read_rds(here::here("data/data_for_analysis.rds")) %>%
  mutate(across(c(epic_cvd, case_diab_caco, fasting, educc3, smk_4cat, alccat,
                  antihyp, lipidlower, ass), as_factor)) %>%
  left_join(new_scores) %>%
  left_join(df_lips)



# Combine phenotype data and scores to create outcome-specific datasets
diab_df <- pheno_df %>%
  filter(ia_subcohort == 1 | case_diab_caco == 1,
         start_time_diab <= stop_time_diab,
         !is.na(start_time_diab)) %>%
  mutate(case_diab_caco =  if_else(case_diab_caco == 1, 1, 0)) %>%
  rename(start = start_time_diab,
         stop = stop_time_diab,
         outcome = case_diab_caco) %>%
  mutate(across(c(cluster_vector, comb_score), ~scale(.) %>%
                  as.numeric()))

cvd_df <- pheno_df %>%
  filter(ia_subcohort == 1 | epic_cvd == 1,
         start_time_cvd <= stop_time_cvd,
         !is.na(start_time_cvd)) %>%
  mutate(epic_cvd =  if_else(epic_cvd == 1, 1, 0)) %>%
  rename(start = start_time_cvd,
         stop = stop_time_cvd,
         outcome = epic_cvd) %>%
  mutate(across(c(cluster_vector, comb_score), ~scale(.) %>%
                  as.numeric()))


cox_fct <- function(.score) {

  diab <-   coxph(Surv(start, stop, outcome) ~
                    eval(rlang::sym(.score)) + sex + cluster(omics_id) + strata(round(age, 0)) +
                    height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
                    smk_4cat + alccat + educc3 + syst + diast,
                  ties = "efron",
                  data = diab_df) %>%
    broom::tidy(exponentiate = TRUE,
                conf.int = TRUE,
                conf.level = 0.95) %>%
    mutate(outcome = "Type 2 diabetes") %>%
    filter(term == "eval(rlang::sym(.score))")

  cvd <-   coxph(Surv(start, stop, outcome) ~
                   eval(rlang::sym(.score)) + sex + cluster(omics_id) + strata(round(age, 0)) +
                   height + waist + aktiv + fasting + antihyp + lipidlower + ass + gj +
                   smk_4cat + alccat + educc3 + syst + diast + prevalent_diab,
                 ties = "efron",
                 data = cvd_df) %>%
    broom::tidy(exponentiate = TRUE,
                conf.int = TRUE,
                conf.level = 0.95) %>%
    mutate(outcome = "CVD") %>%
    filter(term == "eval(rlang::sym(.score))")

  bind_rows(diab, cvd)
}

vars_to_map <- c(cluster_vector, "comb_score") %>%
  set_names(.)


cluster_res <- map_dfr(vars_to_map, cox_fct,
                       .id = "cluster") %>%
  select(-term)



cluster_score_figure <- cluster_res %>%
  mutate(cluster = if_else(cluster == "comb_score",
                           "MLS",
                           str_replace(cluster, "cluster_", "Cluster ") |>
           str_c(" score")) |>
           fct_rev(),
         panel = if_else(cluster == "MLS", 1, 2)) |>
  ggplot() +
  ggforestplot::geom_stripes(aes(y = cluster)) +
  geom_vline(xintercept = 1,
             lty = 2) +
  labs(x = "HR per SD (95%-CI)",
       y = NULL) +
  theme_light(base_family = "RobotoCondensed-Regular") +
  scale_fill_manual(values = c("black", paletteer::paletteer_d("fishualize::Cirrhilabrus_solorensis",
                                                               direction = -1))) +
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
                      y = cluster,
                      xmin = conf.low,
                      xmax = conf.high,
                      fill = cluster),
                  orientation = "y",
                  shape = 22,
                  show.legend = FALSE) +
  facet_grid(panel~outcome,
             scales = "free_y",
             space = "free_y") +
  scale_x_log10(breaks = c(seq(0.6, 1, 0.1),
                           map_dbl(seq(0.6, 1, 0.1), ~1/.) %>% round(digits = 1)))


library(patchwork)

# old multipanel figure
# figure <- network_figure + cluster_score_figure + patchwork::plot_annotation(tag_levels = "A")  +
#   patchwork::plot_layout(ncol = 1,
#                          heights = c(3,1))


# Set path
path_to <- ("/home/fabian/alle_shortcut/!MEP/Projekte/EPIC-Potsdam/Diabetes/Lipidomics/MultiLipidScore/AIP/")

path <- str_c(path_to, "Extended_Figure_9")

# Save as pdf
ggsave(plot = cluster_score_figure,
       glue::glue("{path}.pdf"),
       device = cairo_pdf,
       width = 20,
       height = 10,
       units = "cm")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)
