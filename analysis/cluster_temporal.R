library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(egg)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("temporal_ordering")
dir.create(wd, showWarnings = FALSE)

# set directories, import files ------------------------------------------------
###############################################################################T

function_clusters <- syn("syn21478380") %>%
  read_csv()
pairwise_lfc <- syn("syn21432183") %>%
  read_csv()
meta <- syn("syn21432975") %>%
  read_csv()

condition_meta <- meta %>%
  distinct(condition, ERKi, Time, DOX)

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

# Standardize data -------------------------------------------------------------
###############################################################################T

standardize_data <- function(df, groups, ...) {
  groups_quo <- enquo(groups)
  vars_quo <- enquos(...)
  df %>%
    group_by(!!groups_quo) %>%
    mutate_at(vars(!!!vars_quo), scale, scale = TRUE, center = TRUE) %>%
    ungroup()
}

lfc_long_z <- lfc_long %>%
  standardize_data(gene_id, log2FoldChange) %>%
  inner_join(condition_meta, by = "condition")

# Plotting functions -----------------------------------------------------------
###############################################################################T

plot_single_surface <- function(d, aesthetic = aes(ERKi, Time, fill = log2FoldChange)) {
  fill_quo <- aesthetic[["fill"]]
  abs_max <- quantile(abs(pull(d, !!fill_quo)), .95)
  d %>%
    mutate_at(vars(Time, ERKi), as.factor) %>%
    ggplot(aesthetic) +
    geom_raster() +
    scale_fill_distiller(
      palette = "RdBu", limits = c(-abs_max, abs_max), oob = scales::squish
    ) +
    labs(
      x = "ERK inhibitor (nM)",
      y = "Time (h)",
      fill = "log2 fold change"
    )
}

plot_single_line_by_erki <- function(d) {
  d %>%
    mutate_at(vars(Time, ERKi), as.factor) %>%
    ggplot(aes(Time, log2FoldChange)) +
    geom_line(aes(color = ERKi, group = ERKi)) +
    scale_color_brewer(palette = "Reds")
}

plot_clusters_line_by_erki <- function(d, cluster) {
  cluster_quo <- enquo(cluster)
  d <- d %>%
    mutate_at(vars(Time, ERKi), as.factor)
  avg <- d %>%
    group_by(!!cluster_quo, ERKi, Time) %>%
    summarize(log2FoldChange = mean(log2FoldChange)) %>%
    ungroup()
  d %>%
    ggplot(aes(Time, log2FoldChange)) +
    geom_line(aes(group = gene_id), color = "#444444") +
    geom_line(data = avg, color = "red") +
    facet_grid(vars(ERKi), vars(!!cluster_quo))
}

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
      data = data %>% arrange(!!order_id_quo), size = .8, alpha = 0.1, color = "#CCCCCC"
    )
  p <- p +
    geom_path() +
    geom_point() +
    theme_minimal() +
    scale_color_viridis_c() +
    theme(legend.position = "bottom")
  if (!missing(facet_x) || !missing(facet_y))
    p <- p +
    facet_grid(vars(!!facet_y_quo), vars(!!facet_x_quo), scales = "free")
  p
}


# Temporal clustering k-medoids and hierarchical -------------------------------
###############################################################################T

temporal_distance_per_cluster <- lfc_long_z %>%
  filter(DOX == 1) %>%
  arrange(condition) %>%
  group_nest(gene_id, gene_name) %>%
  inner_join(function_clusters, by = c("gene_id", "gene_name")) %>%
  group_nest(class, direction, class_combined) %>%
  mutate(
    data = map(
      data,
      function(d) {
        n <- nrow(d)
        mat <- matrix(NA_real_, n, n, dimnames = list(d$gene_id, d$gene_id))
        for (i in 1:n) {
          for (j in i:n) {
            mat[[j, i]] <- sqrt(sum((d[["data"]][[i]][["log2FoldChange"]] - d[["data"]][[j]][["log2FoldChange"]])**2))
          }
        }
        as.dist(mat, diag = TRUE)
      }
    )
  )

temporal_clusters_hclust <- temporal_distance_per_cluster %>%
  mutate(
    data = map(
      data,
      hclust,
      method = "average"
    )
  )

plot_single_surface(
  lfc_long %>%
    inner_join(condition_meta, by = "condition") %>%
    filter(gene_name == "JUNB")
)

lfc_long_z %>%
  inner_join(function_clusters, by = c("gene_id", "gene_name")) %>%
  plot_clusters_line_by_erki(class_combined)

temporal_clusters_kmedoid <- temporal_distance_per_cluster %>%
  crossing(k = 2:4) %>%
  mutate(
    cluster_res = map2(
      data, k,
      cluster::pam
    ),
    cluster_df = map(
      cluster_res,
      ~.x %>%
        chuck("clustering") %>%
        enframe("gene_id", "cluster") %>%
        mutate_at(vars(cluster), as.character)
    )
  )

temporal_clusters_kmedoid_surface_plots <- temporal_clusters_kmedoid %>%
  mutate(
    data = map(
      cluster_df,
      function(df) {
        df_trans <- lfc_long_z %>%
          inner_join(df, by = "gene_id") %>%
          group_by(condition, ERKi, Time, DOX, cluster) %>%
          summarize(
            log2FoldChange_mean = mean(log2FoldChange),
            log2FoldChange_sd = sqrt(var(log2FoldChange))
          ) %>%
          ungroup()
        mean_plot <- df_trans %>%
          plot_single_surface(aes(ERKi, Time, fill = log2FoldChange_mean)) +
            facet_wrap(vars(cluster))
        sd_plot <- df_trans %>%
          plot_single_surface(aes(ERKi, Time, fill = log2FoldChange_sd)) +
          facet_wrap(vars(cluster)) +
          labs(fill = "standard deviation")
        ggarrange(mean_plot, sd_plot, draw = FALSE)
      }
    )
  )

wd_clustering <- file.path(wd, "clustering")
dir.create(wd_clustering, showWarnings = FALSE)

write_rds(
  temporal_clusters_hclust,
  file.path(wd_clustering, "temporal_euclidian_hclust.rds")
)

write_rds(
  temporal_clusters_kmedoid,
  file.path(wd_clustering, "temporal_kmedoids.rds")
)

dir.create(file.path(wd_clustering, "kmedoid_plots"), showWarnings = FALSE)

pwalk(
  temporal_clusters_kmedoid_surface_plots,
  function(class_combined, k, data, ...) {
    cowplot::ggsave2(
      file.path(wd_clustering, "kmedoid_plots", paste0("mean_sd_heatmap_", class_combined, "_", k, ".pdf")),
      data
    )
  }
)

# Temporal clustering based on reaching certain relative threshold -------------
###############################################################################T


temporal_ordering_avg_col_stats <- lfc_long %>%
  inner_join(condition_meta, by = "condition") %>%
  filter(DOX == 1) %>%
  inner_join(function_clusters, by = c("gene_id", "gene_name")) %>%
  arrange(gene_id, gene_name, Time, ERKi) %>%
  group_by(gene_id, gene_name, ERKi) %>%
  mutate(
    log2FoldChange_rescaled = log2FoldChange %>%
      {if (.[order(abs(.), decreasing = TRUE)[1]] > 0) scales::rescale(.) else (scales::rescale(-.))},
    log2FoldChange_var = var(log2FoldChange)
  ) %>%
  ungroup() %>%
  group_by(gene_id, gene_name, Time) %>%
  mutate(log2FoldChange_rescaled_mean = weighted.mean(log2FoldChange_rescaled, log2FoldChange_var)) %>%
  ungroup()

temporal_ordering_avg_col <- temporal_ordering_avg_col_stats %>%
  group_by(gene_id, gene_name) %>%
  summarize(
    max_induction = Time[order(log2FoldChange_rescaled_mean, decreasing = TRUE)[1]],
    # In extremely rare cases log2FoldChange_rescaled_mean never reaches 0.5 in
    # any condition, just taking the time point with max induction then
    mid_induction = if (!any(log2FoldChange_rescaled_mean > 0.5)) max_induction else Time[min(which(log2FoldChange_rescaled_mean > 0.5))]
  ) %>%
  ungroup()

plot_single_surface(
  lfc_long %>%
    inner_join(condition_meta, by = "condition") %>%
    filter(gene_name == "FOSB")
)

plot_single_surface(
  temporal_ordering_avg_col_stats %>%
    filter(gene_name == "FOSB") %>%
    mutate(log2FoldChange_rescaled_var = log2FoldChange_rescaled*log2FoldChange_var),
  aes(ERKi, Time, fill = log2FoldChange_rescaled_var)
)

hm_mid_vs_max <- temporal_ordering_avg_col %>%
  count(mid_induction, max_induction) %>%
  mutate_at(vars(mid_induction, max_induction), . %>% as.character() %>%  fct_inseq()) %>%
  ggplot(aes(mid_induction, max_induction, fill = n)) +
    geom_raster()
cowplot::ggsave2(
  file.path(wd, "temporal_ordering_heatmap_mid_vs_max.pdf"),
  hm_mid_vs_max, width = 7, height = 5
)

barplot_mid_vs_max <- temporal_ordering_avg_col %>%
  mutate_at(vars(mid_induction, max_induction), . %>% as.character() %>%  fct_inseq()) %>%
  gather("key", "value", -gene_id, -gene_name) %>%
  mutate_at(vars(value), . %>% as.character() %>%  fct_inseq()) %>%
  ggplot(aes(value, fill = key)) +
    geom_bar(position = "dodge")
cowplot::ggsave2(
  file.path(wd, "temporal_ordering_barplot_mid_vs_max.pdf"),
  barplot_mid_vs_max, width = 5, height = 3
)

write_rds(
  temporal_ordering_avg_col_stats,
  file.path(wd, "temporal_ordering_avg_col_stats.rds")
)

write_csv(
  temporal_ordering_avg_col,
  file.path(wd, "temporal_ordering_avg_col.csv")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

activity <- Activity(
  "Classify genes by the timing of their induction",
  used = c(
    "syn21478380",
    "syn21432183",
    "syn21432975"
  ),
  executed = "https://github.com/clemenshug/erk_senescence/blob/master/analysis/cluster_temporal.R"
)

temporal_ordering_syn <- Folder("temporal_ordering", "syn21432134") %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(wd, "temporal_ordering_avg_col_stats.rds"),
  file.path(wd, "temporal_ordering_avg_col.csv"),
  file.path(wd_clustering, "temporal_euclidian_hclust.rds"),
  file.path(wd_clustering, "temporal_kmedoids.rds")
) %>%
  synStoreMany(temporal_ordering_syn, activity = activity)
