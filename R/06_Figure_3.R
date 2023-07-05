# Functions and packages
source(here::here("R/00_functions.R"))
library(tidyverse)
library(corrr)

# data
data_for_analysis <- read_rds(here::here("data/data_for_analysis.rds")) %>%
  filter(ia_subcohort == 1) %>%
  select(contains("group"), mufa_score, pufa_score, comb_score, sex, age, bmi, waist, height,
         contains("corr"), non_hdl, syst, diast, aktiv) %>%
  rename_food_groups()  %>%
  mutate(across(c(mufa_score, pufa_score, comb_score), ~scale(.) %>%
                  as.numeric()))


# Raw distribution figure -------------------------------------------------

raw_dist <- data_for_analysis %>%
  select(comb_score) %>%
  pivot_longer(everything()) %>%
  ggplot(aes(value)) +
  ggdist::stat_halfeye(adjust = 0.5,
                       justification = -0.07,
                       aes(fill = name),
                       show.legend = FALSE) +
  ggdist::stat_dots(side = "bottom",
                    justification = 1.07,
                    aes(color = name),
                    dotsize = 1.5,
                    show.legend = FALSE) +
  labs(y = NULL,
       x = "Multi-lipid score [Z-scores]") +
  theme_light(base_family = "RobotoCondensed-Regular") +
  scale_x_continuous(breaks = seq(-5, 5, 2.5),
                     limits = c(-5, 5)) +
  scale_color_manual(values = c("MUFA score" = "#ffba00",
                                "Mixed UFA score" = "#00e8ab",
                                "UFA score" = "grey80")) +
  scale_fill_manual(values = c("MUFA score" = "#ffba00",
                               "Mixed UFA score" = "#00e8ab",
                               "UFA score" = "grey80")) +
  theme(strip.background = element_rect(fill = "#333333"),
        strip.text = element_text(color = "white"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = "#333333"))



# Distributions by sex ----------------------------------------------------


fig_sex_diff <- data_for_analysis %>%
  select(comb_score, sex) %>%
  pivot_longer(-sex) %>%
  mutate(sex = if_else(sex == 1, "Men", "Women")) %>%
  ggplot(aes(x = value,
             y = as.factor(sex),
             fill = as.factor(sex))) +
  geom_violin(lwd = 0.2,
              show.legend = FALSE) +
  geom_boxplot(width = 0.2,
               fill = "white",
               outlier.size = 0.8,
               lwd = 0.4,
               color = "#333333",
               outlier.color = "#333333") +
  scale_fill_manual(values = c("Women" = "#e39dff",
                               "Men" = "#3fbea7")) +
  theme_light(base_family = "RobotoCondensed-Regular") +
  labs(x = "Multi-lipid score [Z-scores]",
       y = NULL) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1),
        panel.grid = element_line(color = "grey70"),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "#333333"),
        axis.text = element_text(color = "#333333"))




# Figure on correlation with covars ---------------------------------------

covar_corr <- data_for_analysis %>%
  select(comb_score, bmi, waist, age, non_hdl, corr_chol,
         corr_trigly, corr_hdl, syst, diast) %>%
  correlate(method = "spearman", use = "pairwise.complete.obs", quiet = TRUE) %>%
  focus(comb_score) %>%
  pivot_longer(-term) %>%
  mutate(term = str_replace(term, "bmi", "BMI") %>%
           str_replace("waist", "WC") %>%
           str_replace("syst", "Systolic BP") %>%
           str_replace("diast", "Diastolic BP") %>%
           str_replace("age", "Age") %>%
           str_replace("non_hdl", "Non-HDL-C") %>%
           str_replace("corr_chol", "TC") %>%
           str_replace("corr_trigly", "TG") %>%
           str_replace("corr_hdl", "HDL-C")) %>%
  ggplot(aes(x = value,
             y = fct_reorder(term, value))) +
  scale_y_discrete(position = "right", expand = c(0, 0)) +
  scale_x_continuous(limits = c(-0.6, 0.1), expand = c(0, 0)) +
  geom_vline(xintercept = 0) +
  geom_linerange(aes(xmin = 0, xmax = value),
                 lty = 3,
                 position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5),
             size = 3) +
  labs(x = "Spearman correlation coefficient",
       y = NULL,
       fill = NULL) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1),
        panel.grid = element_line(color = "grey70"),
        panel.grid.minor = element_blank(),
        legend.position = c(0.2, 0.98),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_text(color = "#333333"),
        text = element_text(family = "RobotoCondensed-Regular")) +
  coord_cartesian(clip = "off")  +
  geom_text(data = function(x) filter(x, value < 0),
            aes(x = value - 0.01,
                label = term),
            hjust = 1,
            vjust = 0.5,
            size = 3) +
  geom_text(data = function(x) filter(x, value > 0),
            aes(x = value + 0.01,
                label = term),
            hjust = 0,
            vjust = 0.5,
            size = 3)




# Food correlations -------------------------------------------------------

corr_df <- read_rds(here::here("data/data_for_analysis.rds")) %>%
  filter(ia_subcohort == 1) %>%
  select(contains("group"), comb_score) %>%
  rename_food_groups() %>%
  corrr::correlate(.,method = "spearman",
                   use = "pairwise.complete.obs",
                   quiet = TRUE) %>%
  corrr::focus(comb_score) %>%
  pivot_longer(-term) %>%
  mutate(term = str_to_sentence(term) %>%
           str_replace_all("_", " "))


foods_to_remove <- c("Beer", "Water", "Fruit juice", "Spirits", "Low energy soft drinks", "Coffee",
                     "De caf coffee", "Wine", "Coffee", "Tea")





food_corr <- corr_df %>%
  mutate(highlight = case_when(term == "Butter" ~ "red",
                               term == "Margarine" ~ "blue",
                               TRUE ~ "black")) |>
  # filter(!term %in% foods_to_remove) %>%
  # group_by(value < 0) %>%
  # slice_max(abs(value),
  #           n = 5) %>%
  ggplot(aes(y = value,
             x = fct_reorder(term, value),
             color = highlight)) +
  scale_x_discrete(expand = c(0, 0.5)) +
  scale_y_continuous(limits = c(-0.26, 0.25), expand = c(0, 0)) +
  geom_hline(yintercept = 0) +
  geom_linerange(aes(ymin = 0, ymax = value),
                 lty = 3,
                 position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5),
             size = 3) +
  labs(y = "Spearman correlation coefficient",
       x = NULL,
       fill = NULL) +
  scale_color_identity() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(size = 0.1),
        panel.grid = element_line(color = "grey70"),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(color = "#333333"),
        axis.text.x = element_blank(),
        text = element_text(family = "RobotoCondensed-Regular")) +
  coord_cartesian(clip = "off",
                  expand = TRUE) +
  geom_text(data = function(x) slice_max(x, value),
            aes(x = term,
                label = term),
            hjust = 1,
            nudge_x = -1,
            vjust = 0.5,
            size = 5) +
  geom_text(data = function(x) slice_min(x, value),
            aes(x = term,
                label = term),
            hjust = 0,
            nudge_x = 1,
            vjust = 0.5,
            size = 5)



# Combine figures into panel figure
library(patchwork)


layout <- "
AABB
CCDD
"


figure <- raw_dist + fig_sex_diff + covar_corr + food_corr +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = "A")


# Set path
path <- here::here("doc", "Figures", "Figure_3")

# Save as pdf
ggsave(plot = figure,
       glue::glue("{path}.pdf"),
       device = cairo_pdf,
       width = 22,
       height = 18,
       units = "cm")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)

