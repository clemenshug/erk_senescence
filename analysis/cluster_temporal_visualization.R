library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(pheatmap)
library(RColorBrewer)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("temporal_ordering")
dir.create(wd, showWarnings = FALSE)

# set directories, import files ------------------------------------------------
###############################################################################T

surface_fit <- syn("syn21444486") %>%
  read_csv()

temporal_ordering <- syn("syn21536903") %>%
  read_csv()

temporal_ordering_stats <- syn("syn21537190") %>%
  read_csv()

meta <- syn("syn21432975") %>%
  read_csv()

function_clusters <- syn("syn21478380") %>%
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
    directed,
    Time,
    log2FoldChange_rescaled_mean
  ) %>%
  inner_join(temporal_ordering, by = c("gene_id", "gene_name", "directed"))

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

plot_induction_heatmap <- function(
  df, row_id, row_order, col_id, col_id_order, fill, ...
) {
# plot_induction_heatmap <- function(df, row_id, col_id, fill) {
  row_id_quo <- enquo(row_id)
  row_order_quo <- enquo(row_order)
  col_id_quo <- enquo(col_id)
  col_order_quo <- enquo(col_id_order)
  fill_quo <- enquo(fill)
  row_id_order <- df %>%
    arrange(!!!eval(rlang::get_expr(row_order_quo), envir = df)) %>%
    pull(!!row_id_quo) %>%
    unique()
  col_id_order <- df %>%
    arrange(!!!eval(rlang::get_expr(col_order_quo), envir = df)) %>%
    pull(!!col_id_quo) %>%
    unique()
  df <- df %>%
    mutate(
      !!row_id_quo := factor(!!row_id_quo, levels = row_id_order),
      !!col_id_quo := factor(!!col_id_quo, levels = col_id_order)
    )
  abs_max <- df %>%
    pull(!!fill_quo) %>%
    abs() %>%
    quantile(.90, names = FALSE, na.rm = TRUE) %>%
    round(1)
  if (abs_max < 0.01) abs_max <- 0.01
  breaks <- seq(-abs_max, abs_max, by = signif(abs_max/10, 1))
  cmap <- rev(colorRampPalette(brewer.pal(7, "RdBu"))(length(breaks) - 1))
  # message(breaks)
  mat <- df %>%
    select(!!row_id_quo, !!col_id_quo, !!fill_quo) %>%
    spread(!!col_id_quo, !!fill_quo) %>%
    column_to_rownames(rlang::as_name(row_id_quo)) %>%
    as.matrix()
  # browser()
  pheatmap(
    mat,
    color = cmap,
    breaks = breaks,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    draw = FALSE,
    ...
  )
}



temporal_ordering_heatmap_data <- temporal_ordering_stats %>%
  mutate_at(vars(ERKi), as.character) %>%
  inner_join(temporal_ordering, by = c("gene_id", "gene_name", "ERKi", "directed")) %>%
  distinct(
    gene_id,
    directed,
    direction_max = direction_max.x,
    log2FoldChange,
    log2FoldChange_norm = log2FoldChange_rescaled * log2FoldChange_var_erki * if_else(direction_max.x == "pos", 1, -1),
    ERKi,
    Time,
    mid_induction,
    max_induction,
    variance
  )

temporal_ordering_heatmap_row_ann <- temporal_ordering %>%
  select(gene_id, ERKi, directed, max_induction, mid_induction) %>%
  mutate_at(
    vars(ends_with("induction")),
    ~factor(as.character(.x), levels = as.character(sort(unique(.x))))
  ) %>%
  inner_join(select(function_clusters, gene_id, class_combined), by = c("gene_id")) %>%
  group_nest(ERKi, directed, .key = "ann") %>%
  mutate(
    ann = map(
      ann,
      ~.x %>%
        column_to_rownames("gene_id") %>%
        as.data.frame()
    )
  )


temporal_ordering_heatmaps <- temporal_ordering_heatmap_data %>%
  group_nest(ERKi, directed) %>%
  inner_join(temporal_ordering_heatmap_row_ann, by = c("ERKi", "directed")) %>%
  mutate(
    data = map2(
      data, ann,
      ~plot_induction_heatmap(
        .x, gene_id, list(direction_max, mid_induction, desc(variance)), Time, list(Time), log2FoldChange,
        show_rownames = FALSE, annotation_row = .y
      )
    )
  )

pwalk(
  temporal_ordering_heatmaps,
  function(ERKi, directed, data, ...) {
    ggsave(
      file.path(wd, paste0("temporal_ordering_heatmap_erki_", ERKi, "_", directed, ".pdf")),
      data,
      width = 4, height = 5
    )
  }
)

# Scatter plot of gene ranks at different ERK concentrations -------------------
###############################################################################T

scatter_mid_max <- temporal_ordering %>%
  arrange(ERKi, mid_induction, max_induction, desc(variance)) %>%
  filter(ERKi %in% c("0", "1000")) %>%
  mutate(ERKi = paste0("rank_ERKi_", ERKi)) %>%
  group_by(ERKi, directed) %>%
  mutate(rank = 1:n()) %>%
  ungroup() %>%
  group_by(gene_name, directed) %>%
  summarize(
    rank = list(set_names(rank, ERKi)),
    variance = prod(variance)
  ) %>%
  ungroup() %>%
  inner_join(function_clusters, by = "gene_name") %>%
  unnest_wider(rank) %>%
  ggplot(aes(rank_ERKi_0, rank_ERKi_1000, alpha = variance, color = class_combined)) +
    geom_point() +
    scale_alpha_continuous(trans = "log10") +
    facet_grid(vars(directed), vars(class_combined))

ggsave(
  file.path(wd, "temporal_ordering_rank_erki_0_vs_erki_1000.pdf"),
  scatter_mid_max, width = 30, height = 8
)

