library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(pheatmap)
library(RColorBrewer)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("clustering")
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

dose_clusters <- syn("syn21576614") %>%
  read_csv()

condition_meta <- meta %>%
  distinct(condition, ERKi, Time, DOX)

pairwise_lfc <- syn("syn21432183") %>%
  read_csv()

pairwise_padj <- syn("syn21432184") %>%
  read_csv()

padj_long <- pairwise_padj  %>%
  mutate_at(vars(starts_with("time")), replace_na, 0) %>%
  gather("condition", "padj", -gene_name, -gene_id)

deseq_lfc <- syn("syn21432183") %>%
  read_csv() %>%
  gather("condition", "log2_FC", -starts_with("gene")) %>%
  drop_na()


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

pca_values <- syn("syn21444456") %>%
  read_csv() %>%
  gather("PC", "PC_value", -condition)

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

deseq_lfc_scaled <- deseq_lfc %>%
  group_by(gene_id) %>%
  mutate_at(vars(log2_FC), ~scale(.x, center = FALSE)[, 1]) %>%
  ungroup()


# Use different color scheme
# So that there's no confusion

consensus_trajectories_plot <- plot_cluster_trajectories(
  consensus_trajectories_data %>%
    filter(!consensus %in% c("none", "no_response_0")) %>%
    mutate(consensus = cluster_names[consensus]),
  PC_value, log2_FC, condition, PC1_order, gene_id, PC, consensus,
  all_traces = TRUE, rescale_limits = 0.4
) +
  guides(color = FALSE)


cluster_classes_long <- dose_clusters %>%
  select(-gene_name, -consensus_n) %>%
  pivot_longer(
    c(function_fitting, k_medoids, linear_model, consensus),
    names_to = "method",
    values_to = "dose_class"
  ) %>%
  filter(!dose_class %in% c("no_response_0", "none")) %>%
  mutate(
    direction = str_extract(dose_class, "[+-]"),
    dose_class = str_replace(dose_class, "_[+-]", "")
  )

cluster_traces_data <- condition_meta %>%
  filter(Time == "24", condition != "time0dox0conc0") %>%
  inner_join(
    pca_values %>%
      filter(PC %in% c("PC1", "PC2")),
    by = "condition"
  ) %>%
  crossing(
    cluster_classes_long
  ) %>%
  inner_join(
    deseq_lfc,
    by = c("gene_id", "condition")
  ) %>%
  inner_join(
    pca_values %>%
      filter(PC == "PC1") %>%
      dplyr::select(condition, PC1_order = PC_value),
    by = "condition"
  ) %>%
  mutate(
    class_combined = paste(dose_class, direction, sep = "_")
  )


cluster_traces <- cluster_traces_data %>%
  group_nest(method) %>%
  rowwise() %>%
  mutate(
    plot = plot_cluster_trajectories(
      group_by(data, gene_id) %>%
        mutate_at(vars(log2_FC), scale, scale = TRUE, center = TRUE) %>%
        ungroup(),
      PC_value, log2_FC, condition, PC1_order, gene_id, PC, class_combined,
      all_traces = TRUE
    ) %>%
      list()
  ) %>%
  ungroup()

pwalk(
  cluster_traces,
  function(method, plot, ...) {
    cowplot::ggsave2(
      file.path(wd, paste0("cluster_traces_", method, ".pdf")),
      plot,
      width = 5, height = 10
    )
  }
)
