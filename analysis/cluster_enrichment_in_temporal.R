library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(RColorBrewer)
library(ggforce)
library(viridis)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("umap")
dir.create(wd, recursive = TRUE, showWarnings = FALSE)

# set directories, import files ------------------------------------------------
###############################################################################T

meta <- syn("syn21432975") %>%
  read_csv()

condition_meta <- meta %>%
  distinct(condition, DMSO, ERKi,  Time, DOX)

cluster_overlap <- syn("syn21576614") %>%
  read_csv()

cluster_names <- syn("syn21567677") %>%
  read_csv() %>%
  bind_rows(list(class_combined = "unknown", class_name = "Unknown")) %>%
  mutate_at(vars(class_name), . %>% as.factor() %>% fct_inorder()) %>%
  with(set_names(class_name, class_combined))

temporal_ordering <- syn("syn21536903") %>%
  read_csv()

# Calculate enrichments --------------------------------------------------------
###############################################################################T

calc_enrichment <- function(x, y, t) {
  v <- t[y, x]
  x_o <- sum(t[y, !colnames(t) %in% x])
  y_o <- sum(t[!rownames(t) %in% y, x])
  xy_o <- sum(t[!rownames(t) %in% y, !colnames(t) %in% x])
  m <- matrix(c(v, x_o, y_o, xy_o), ncol = 2, byrow = TRUE)
  # browser()
  fisher.test(m)
}

enrichment <- temporal_ordering %>%
  filter(directed == "directed", ERKi %in% c("0", "1000")) %>%
  inner_join(
    cluster_overlap,
    by = c("gene_id", "gene_name")
  ) %>%
  group_by(ERKi) %>%
  summarize(
    cont_table = table(mid_induction, consensus) %>%
      list()
  ) %>%
  ungroup() %>%
  mutate(
    res = map(
      cont_table,
      function(t) {
        crossing(
          consensus = colnames(t),
          mid_induction = rownames(t)
        ) %>%
          mutate(
            fisher_res = map2(
              consensus, mid_induction,
              calc_enrichment,
              t
            ) %>%
              map(broom::tidy)
          ) %>%
          unnest(fisher_res) %>%
          mutate(padj = p.adjust(p.value, method = "fdr"))
      }
    )
  )

enrichment_heatmap <- enrichment %>%
  unnest(res) %>%
  mutate(
    consensus = cluster_names[consensus],
    mid_induction = as.integer(mid_induction) %>%
      factor(levels = sort(unique(.)))
  ) %>%
  ggplot(
    aes(
      mid_induction,
      consensus,
      fill = -log2(estimate),
      label = cut(padj, breaks = c(-Inf, 0.001, 0.05, Inf), labels = c("**", "*", ""))
    )
  ) +
  geom_tile() +
  geom_text() +
  scale_fill_distiller(palette = "RdBu", limits = c(-2.5, 2.5), oob = scales::squish) +
  facet_wrap(vars(ERKi))

enrichment_barplot <- enrichment %>%
  select(-res) %>%
  mutate(cont_table = map(cont_table, ~.x %>% as_tibble())) %>%
  unnest(cont_table) %>%
  group_by(consensus, ERKi) %>%
  mutate(percent = 100*n/sum(n)) %>%
  ungroup() %>%
  mutate(
    consensus = cluster_names[consensus] %>%
      fct_recode(NULL = "other"),
    mid_induction = as.integer(mid_induction) %>%
      factor(levels = sort(unique(.)))
  ) %>%
  drop_na(consensus) %>%
  ggplot(aes(consensus, percent, fill = mid_induction)) +
    geom_col(position = "dodge") +
    scale_fill_brewer(palette = "Reds") +
    facet_wrap(vars(ERKi)) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

enrichment_barplot <- enrichment %>%
  select(-res) %>%
  mutate(cont_table = map(cont_table, ~.x %>% as_tibble())) %>%
  unnest(cont_table) %>%
  group_by(consensus, ERKi) %>%
  mutate(percent = 100*n/sum(n)) %>%
  ungroup() %>%
  mutate(
    consensus = cluster_names[consensus] %>%
      fct_recode(NULL = "other"),
    mid_induction = as.integer(mid_induction) %>%
      factor(levels = sort(unique(.)))
  ) %>%
  drop_na(consensus) %>%
  ggplot(aes(mid_induction, percent, fill = ERKi)) +
  geom_col(position = "dodge") +
  # geom_line(aes(group = ERKi)) +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(vars(consensus))
  # theme(axis.text.x = element_text(angle = 30, hjust = 1))

