library(tidyverse)
library(tidygraph)
library(ggraph)


matching_table <- read_rds(here::here("doc/matching_for_predimed.rds")) %>%
  mutate(across(c(lipid_display, within_class_fa_sum), ~str_to_lower(.) %>%
                  str_replace_all("dg", "dag") %>%
                  str_replace_all("tg", "tag") %>%
                  str_replace_all("\\(", "_") %>%
                  str_remove_all("\\)") %>%
                  str_replace_all(":", "_") %>%
                  str_replace_all("pe_", "peonly_") %>%
                  str_replace_all("pep_", "pe_p_") %>%
                  str_replace_all("peo_", "pe_o_")))

matching_table %>%
  count(within_class_fa_sum)


no_match_metab <- matching_table %>%
  group_by(within_class_fa_sum) %>%
  summarize(matches = paste(metabolite, collapse = ",")) %>%
  mutate(temp = parse_number(matches)) %>%
  filter(is.na(temp))


root_to_weights <- matching_table %>%
  select(to = within_class_fa_sum) %>%
  mutate(from = "MLS")

temp <- matching_table %>%
  na.omit() %>%
  select(from = within_class_fa_sum,
         to = lipid_display) %>%
  bind_rows(root_to_weights, .) %>%
  as_tbl_graph() %>%
  mutate(miss = if_else(name %in% no_match_metab$within_class_fa_sum,
                        "no match",
                        "with match"),
         lips = str_to_upper(name) %>%
           str_replace("TAG", "TG") %>%
           str_replace("MAG", "MG") %>%
           str_replace("DAG", "DG") %>%
           str_replace("LACCER", "LacCer") %>%
           str_replace("HEXCER", "HexCer") %>%
           str_replace("DHCER", "dhCer") %>%
           str_replace("^CER", "Cer") %>%
           str_remove("ONLY") %>%
           str_replace("PE_O", "PEO") %>%
           str_replace("PE_P", "PEP"))

temp2 <- temp %>%
  activate(nodes) %>%
  as_tibble() %>%
  separate(lips, remove = FALSE, into = c("base", "l", "db"), fill = "right") %>%
  mutate(display_name = if_else(base == "MLS", "MLS", str_c(base, "(", l, ":", db, ")")))


figure <- temp %>%
  activate(nodes) %>%
  left_join(temp2) %>%
  ggraph("tree") +
  geom_edge_diagonal() +
  labs(fill = NULL) +
  geom_node_label(aes(label = display_name),
                  size = 4) +
  geom_node_label(aes(fill = miss,
                      label = display_name),
                  size = 4,
                  data = function(x) filter(x, name %in% root_to_weights$to)) +
  scale_fill_manual(values = c("red", "steelblue")) +
  theme_void() +
  coord_flip() +
  scale_y_reverse() +
  scale_x_reverse() +
  theme(legend.position = c(0.9, 0.1),
        legend.title = element_blank(),
        legend.text = element_text(size = 15)) +
  guides(fill = guide_legend(order = 2, override.aes = aes(label = "",
                                                           size = 15)),
         shape = "none") +
  # annotate("text", x = -2, y = 2, label = "MLS") +
  annotate("text", x = -41, y = 1, label = "Lipids in score", size = 8) +
  annotate("text", x = -40, y = 0.3, label = "Predictors of missing\n score lipids", size = 8)



# Set path
path_to <- ("/home/fabian/alle_shortcut/!MEP/Projekte/EPIC-Potsdam/Diabetes/Lipidomics/MultiLipidScore/AIP/")

path <- str_c(path_to, "Extended_Figure_4")

# Save as pdf
ggsave(plot = figure,
       glue::glue("{path}.pdf"),
       device = cairo_pdf,
       height = 40,
       width = 25,
       units = "cm",
       bg = "white")

# Convert to png and save
pdftools::pdf_convert(pdf = glue::glue("{path}.pdf"),
                      filenames = glue::glue("{path}.png"),
                      format = "png",
                      dpi = 400)
