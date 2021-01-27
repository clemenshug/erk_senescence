library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(RColorBrewer)
library(ggforce)
library(viridis)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("clustering", "overlap")
dir.create(wd, recursive = TRUE, showWarnings = FALSE)

# set directories, import files ------------------------------------------------
###############################################################################T

deseq_lfc <- syn("syn21432183") %>%
  read_csv() %>%
  gather("condition", "log2_FC", -starts_with("gene")) %>%
  drop_na()

meta <- syn("syn21432975") %>%
  read_csv()

condition_meta <- meta %>%
  distinct(condition, DMSO, ERKi,  Time, DOX)

cluster_overlap <- syn("syn21576614") %>%
  read_csv()

pca_values <- syn("syn21444456") %>%
  read_csv() %>%
  gather("PC", "PC_value", -condition)

cluster_names <- syn("syn21567677") %>%
  read_csv() %>%
  bind_rows(list(class_combined = "unknown", class_name = "Unknown")) %>%
  mutate_at(vars(class_name), . %>% as.factor() %>% fct_inorder()) %>%
  with(set_names(class_name, class_combined))


# Plot cluster trajectories ----------------------------------------------------
###############################################################################T


plot_cluster_trajectories <- function(
  data, x, y, condition_id, order_id, trajectory_id, facet_x = NULL, facet_y = NULL,
  all_traces = FALSE, ci = NULL, rescale_limits = NULL
) {
  x_quo <- enquo(x)
  y_quo <- enquo(y)
  condition_id_quo <- enquo(condition_id)
  trajectory_id_quo <- enquo(trajectory_id)
  order_id_quo <- enquo(order_id)
  facet_x_quo <- enquo(facet_x)
  facet_y_quo <- enquo(facet_y)
  # browser()
  averages <- data %>%
    group_by(!!order_id_quo, !!condition_id_quo, !!facet_x_quo, !!facet_y_quo) %>%
    summarize_at(vars(!!x_quo, !!y_quo), mean) %>%
    ungroup() %>%
    arrange(!!order_id_quo)
  p <- ggplot(averages, aes(!!x_quo, !!y_quo, color = !!order_id_quo))
  if (all_traces) {
    p <- p +
      geom_path(
        aes(group = !!trajectory_id_quo),
        data = arrange(data, !!order_id_quo), size = .8, alpha = 0.1, color = "#CCCCCC"
      )
    if (!is.null(rescale_limits))
      p <- p +
        scale_y_continuous(limits = function(x) x*rescale_limits, oob = function(x, y) x)
  }
  if (!is.null(ci)) {
    intervals <- data %>%
      group_by(!!order_id_quo, !!condition_id_quo, !!facet_x_quo, !!facet_y_quo) %>%
      summarize(
        x = mean(!!x_quo),
        ymin = quantile(!!y_quo, 1 - ci),
        ymax = quantile(!!y_quo, ci)
      ) %>%
      ungroup() %>%
      arrange(!!order_id_quo)
    p <- p +
      geom_ribbon(
        aes(x = x, y = NULL, ymin = ymin, ymax = ymax, color = NULL),
        data = intervals,
        orientation = "x"
      )
  }
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

deseq_lfc_scaled <- deseq_lfc %>%
  group_by(gene_id) %>%
  mutate_at(vars(log2_FC), ~scale(.x, center = FALSE)[, 1]) %>%
  ungroup()

consensus_trajectories_data <- condition_meta %>%
  filter(Time == "24", condition != "time0dox0conc0") %>%
  inner_join(
    pca_values %>%
      filter(PC %in% c("PC1", "PC2")),
    by = "condition"
  ) %>%
  crossing(
    cluster_overlap
  ) %>%
  inner_join(
    deseq_lfc,
    by = c("gene_id", "gene_name", "condition")
  ) %>%
  inner_join(
    pca_values %>%
      filter(PC == "PC1") %>%
      dplyr::select(condition, PC1_order = PC_value),
    by = "condition"
  )


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

ggsave(
  file.path(wd, "consensus_clusters_all_traces.pdf"),
  consensus_trajectories_plot,
  width = 4, height = 12
)


consensus_trajectories_plot_average <- plot_cluster_trajectories(
  consensus_trajectories_data %>%
    filter(!consensus %in% c("none", "no_response_0")),
  PC_value, log2_FC, condition, PC1_order, gene_id, PC, consensus,
  all_traces = FALSE
) +
  guides(color = FALSE)

# Use different color scheme
# So that there's no confusion

ggsave(
  file.path(wd, "consensus_clusters_average_traces.pdf"),
  consensus_trajectories_plot_average,
  width = 3, height = 10
)

consensus_trajectories_plot_ci <- plot_cluster_trajectories(
  consensus_trajectories_data %>%
    filter(!consensus %in% c("none", "no_response_0")),
  PC_value, log2_FC, condition, PC1_order, gene_id, PC, consensus,
  all_traces = FALSE, ci = 0.95
) +
  guides(color = FALSE)

# Use different color scheme
# So that there's no confusion

ggsave(
  file.path(wd, "consensus_clusters_average_traces.pdf"),
  consensus_trajectories_plot_average,
  width = 3, height = 10
)

