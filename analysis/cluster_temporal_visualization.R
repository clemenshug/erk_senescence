library(tidyverse)
library(ssh)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("temporal_ordering")
dir.create(wd, showWarnings = FALSE)

# set directories, import files ------------------------------------------------
###############################################################################T

surface_fit <- syn("syn21444486") %>%
  read_csv()

temporal_ordering <- syn("syn21521367") %>%
  read_csv()

temporal_ordering_stats <- syn("syn21521366") %>%
  read_rds()

meta <- syn("syn21432975") %>%
  read_csv()

condition_meta <- meta %>%
  distinct(condition, ERKi, Time, DOX)

pairwise_lfc <- syn("syn21432183") %>%
  read_csv()

lfc_long <- pairwise_lfc %>%
  # Remove genes that are NA in all conditions
  filter(
    !select(., starts_with("time")) %>%
      as.matrix() %>%
      is.na() %>%
      apply(1, all)
  ) %>%
  mutate_at(vars(starts_with("time")), replace_na, 0) %>%
  gather("condition", "log2FoldChange", -gene_name, -gene_id)

# Visualize induction clusters -------------------------------------------------
###############################################################################T

plot_cluster_trajectories <- function(
  data, x, y, order_id, trajectory_id, facet_x = NULL, facet_y = NULL,
  all_traces = FALSE
) {
  x_quo <- enquo(x)
  y_quo <- enquo(y)
  order_id_quo <- enquo(order_id)
  trajectory_id_quo <- enquo(trajectory_id)
  facet_x_quo <- enquo(facet_x)
  facet_y_quo <- enquo(facet_y)
  averages <- data %>%
    group_by(!!x_quo, !!facet_x_quo, !!facet_y_quo) %>%
    summarize_at(vars(!!y_quo), mean) %>%
    ungroup() %>%
    arrange(!!order_id_quo)
  p <- ggplot(averages, aes(!!x_quo, !!y_quo, color = !!order_id_quo))
  if (all_traces)
    p <- p +
    geom_path(
      aes(group = !!trajectory_id_quo),
      data = data %>% arrange(!!order_id_quo), size = .8, alpha = 0.2, color = "#CCCCCC"
    )
  p <- p +
    geom_path() +
    geom_point() +
    theme_minimal() +
    scale_color_viridis_c() +
    theme(legend.position = "bottom")
  if (!missing(facet_x) && !missing(facet_y))
    p <- p +
      facet_grid(vars(!!facet_y_quo), vars(!!facet_x_quo), scales = "free")
  else (!missing(facet_x))
    p <- p +
      facet_wrap(vars(!!facet_x_quo), scales = "free")
  p
}

temporal_ordering_avg <- temporal_ordering_stats %>%
  distinct(
    gene_id,
    gene_name,
    Time,
    log2FoldChange_rescaled_mean
  ) %>%
  inner_join(temporal_ordering, by = c("gene_id", "gene_name"))

temporal_ordering_traces <-  plot_cluster_trajectories(
  temporal_ordering_avg %>%
    filter(ERKi == "all") %>%
    distinct(gene_id, Time, max_induction, mid_induction, log2FoldChange_rescaled_mean) %>%
    gather("induction", "induction_time", mid_induction, max_induction) %>%
    mutate_at(vars(induction_time), . %>% as.character %>% fct_reorder(as.numeric(.))),
  Time, log2FoldChange_rescaled_mean, Time, gene_id,
  facet_x = induction, facet_y = induction_time, all_traces = TRUE
) +
  labs(y = "Expression induction") +
  facet_grid(vars(induction), vars(induction_time))

cowplot::ggsave2(
  file.path(wd, "temporal_ordering_all_erki_traces.pdf"),
  temporal_ordering_traces, width = 8, height = 5
)

# Ordered heatmaps of gene induction -------------------------------------------
###############################################################################T

# lfc_long_clusters <- lfc_long %>%
#   inner_join(condition_meta, by = "condition") %>%
#   mutate_at(vars(ERKi), as.character) %>%
#   inner_join(
#     temporal_ordering %>%
#       # rename(induction_ERKi = ERKi) %>%
#       filter(ERKi %in% c("0", "1000")),
#     by = c("gene_id", "gene_name", "ERKi")
#   )

# plot_induction_heatmap <- function(df, row_id, row_order, col_id, col_id_order, fill) {
plot_induction_heatmap <- function(df, row_id, col_id, fill) {
  row_id_quo <- enquo(row_id)
  # row_order_quo <- map(row_order, enquo)
  col_id_quo <- enquo(col_id)
  # col_order_quo <- enquo(col_id_order)
  fill_quo <- enquo(fill)
  # browser()
  # row_id_order <- df %>%
  #   arrange(!!row_order_quo) %>%
  #   pull(!!row_id_quo) %>%
  #   unique()
  # col_id_order <- df %>%
  #   arrange(!!col_order_quo) %>%
  #   pull(!!col_id_quo) %>%
  #   unique()
  # df <- df %>%
  #   mutate(
  #     !!row_id_quo := factor(!!row_id_quo, levels = row_id_order),
  #     !!col_id_quo := factor(!!col_id_quo, levels = col_id_order)
  #   )
  abs_max <- df %>%
    pull(!!fill_quo) %>%
    abs() %>%
    quantile(.90, names = FALSE, na.rm = TRUE) %>%
    round(1)
  # breaks <- seq(-abs_max, abs_max, by = 0.1)
  # cmap <- rev(colorRampPalette(brewer.pal(7, "RdBu"))(length(breaks)))
  df %>%
    ggplot(aes(!!col_id_quo, !!row_id_quo, fill = !!fill_quo)) +
      geom_raster() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      ) +
      scale_fill_distiller(palette = "RdBu", limits = c(-abs_max, abs_max))
}


hm_1000 <- temporal_ordering_stats %>%
  mutate_at(vars(ERKi), as.character) %>%
  inner_join(temporal_ordering, by = c("gene_id", "gene_name", "ERKi")) %>%
  filter(ERKi %in% c("1000")) %>%
  distinct(gene_id, log2FoldChange, ERKi, Time, mid_induction, max_induction) %>%
  group_by(gene_id, ERKi) %>%
  mutate(
    # pos_neg = filter(., Time == 24) %>%
    #   chuck("log2FoldChange") %>%
    #   {if (. > 0) "pos" else "neg"}
    pos_neg = if (median(log2FoldChange) > 0) "pos" else "neg"
  ) %>%
  ungroup() %>%
  mutate(
    Time = factor(
      Time,
      levels = arrange(., Time) %>%
        pull(Time) %>%
        unique()
    ),
    gene_id = factor(
      gene_id,
      levels = arrange(., pos_neg, desc(mid_induction), desc(max_induction)) %>%
        pull(gene_id) %>%
        unique()
    )
  ) %>%
  plot_induction_heatmap(
    gene_id, Time, log2FoldChange
  ) +
  ggtitle("ERKi 1000nM")

hm_0 <- temporal_ordering_stats %>%
  mutate_at(vars(ERKi), as.character) %>%
  inner_join(temporal_ordering, by = c("gene_id", "gene_name", "ERKi")) %>%
  filter(ERKi %in% c("0")) %>%
  distinct(gene_id, log2FoldChange, ERKi, Time, mid_induction, max_induction) %>%
  group_by(gene_id, ERKi) %>%
  mutate(
    # pos_neg = filter(., Time == 24) %>%
    #   chuck("log2FoldChange") %>%
    #   {if (. > 0) "pos" else "neg"}
    pos_neg = if (median(log2FoldChange) > 0) "pos" else "neg"
  ) %>%
  ungroup() %>%
  mutate(
    Time = factor(
      Time,
      levels = arrange(., Time) %>%
        pull(Time) %>%
        unique()
    ),
    gene_id = factor(
      gene_id,
      levels = arrange(., pos_neg, desc(mid_induction), desc(max_induction)) %>%
        pull(gene_id) %>%
        unique()
    )
  ) %>%
  plot_induction_heatmap(
    gene_id, Time, log2FoldChange
  ) +
  ggtitle("ERKi 0nM")

combined <- egg::ggarrange(
  hm_0,
  hm_1000,
  draw = FALSE, ncol = 2
)
ggsave(file.path(wd, "temporal_ordering_heatmap.pdf"), combined, width = 6, height = 10)
ggsave(file.path(wd, "temporal_ordering_heatmap.png"), combined, width = 6, height = 10)

# Scatter plot of gene ranks at different ERK concentrations -------------------
###############################################################################T

scatter_mid_max <- temporal_ordering %>%
  arrange(mid_induction, max_induction) %>%
  filter(ERKi %in% c("0", "1000")) %>%
  mutate(ERKi = paste0("rank_ERKi_", ERKi)) %>%
  group_by(ERKi) %>%
  mutate(rank = 1:n()) %>%
  ungroup() %>%
  group_by(gene_name) %>%
  summarize(
    rank = list(set_names(rank, ERKi)),
    variance = prod(variance)
  ) %>%
  ungroup() %>%
  unnest_wider(rank) %>%
  ggplot(aes(rank_ERKi_0, rank_ERKi_1000, alpha = variance)) +
    geom_point()

ggsave(
  file.path(wd, "temporal_ordering_rank_erki_0_vs_erki_1000.pdf"),
  scatter_mid_max
)

