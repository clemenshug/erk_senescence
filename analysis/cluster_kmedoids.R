library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(polynom)
library(broom)
library(furrr)
library(cluster)
library(SimilarityMeasures)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("clustering", "kmedoids")
dir.create(wd, recursive = TRUE, showWarnings = FALSE)

syn_clustering <- "syn21432135"

# set directories, import files ------------------------------------------------
###############################################################################T

pairwise_lfc <- syn("syn21432183") %>%
  read_csv()
meta <- syn("syn21432975") %>%
  read_csv()
pca_values <- syn("syn21444456") %>%
  read_csv()
surface_fit <- syn("syn21444486") %>%
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

# Add PC info to log2FoldChange ------------------------------------------------
###############################################################################T

pc_lfc_long <- lfc_long %>%
  filter(gene_id %in% surface_fit$gene_id) %>%
  inner_join(select(pca_values, condition, PC1, PC2), by = "condition")
  # inner_join(select(condition_merta, condition, PC1, PC2, erk_range), by = "condition")


# Plotting funcs ---------------------------------------------------------------
###############################################################################T

plot_cluster_trajectories <- function(
  data, x, y, condition_id, order_id, trajectory_id, facet_x = NULL, facet_y = NULL,
  all_traces = FALSE
) {
  x_quo <- enquo(x)
  y_quo <- enquo(y)
  condition_id_quo <- enquo(condition_id)
  trajectory_id_quo <- enquo(trajectory_id)
  order_id_quo <- enquo(order_id)
  facet_x_quo <- enquo(facet_x)
  facet_y_quo <- enquo(facet_y)
  averages <- data %>%
    group_by(!!order_id_quo, !!condition_id_quo, !!facet_x_quo, !!facet_y_quo) %>%
    summarize_at(vars(!!x_quo, !!y_quo), mean) %>%
    ungroup() %>%
    arrange(!!order_id_quo)
  p <- ggplot(averages, aes(!!x_quo, !!y_quo, color = !!order_id_quo))
  if (all_traces)
    p <- p +
    geom_path(
      aes(group = !!trajectory_id_quo),
      data = arrange(data, !!order_id_quo), size = .8, alpha = 0.1, color = "#CCCCCC"
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

# Distance metrics -------------------------------------------------------------
###############################################################################T

standardize_data <- function(df, groups, ...) {
  groups_quo <- enquo(groups)
  vars_quo <- enquos(...)
  df %>%
    group_by(!!groups_quo) %>%
    mutate_at(vars(!!!vars_quo), scale, scale = TRUE, center = TRUE) %>%
    ungroup()
}

distance_euclidian <- function(df, groups, by, ...) {
  groups_sym <- enquo(groups)
  by_sym <- enquo(by)
  vars_quo  <- enquos(...)
  l <- df %>%
    arrange(!!groups_sym, !!by_sym) %>%
    group_nest(!!groups_sym) %>%
    mutate(
      data = map(
        data,
        ~dplyr::select(.x, !!!vars_quo) %>%
          as.matrix()
      )
    )
  n <- nrow(l)
  dist <- matrix(NA_real_, n, n)
  for (i in seq_len(n)) {
    for (j in seq(i, n)) {
      dist[i, j] <- sum(sqrt(rowSums((l$data[[i]] - l$data[[j]])**2)))
    }
  }
  rownames(dist) <- pull(l, !!groups_sym)
  colnames(dist) <- pull(l, !!groups_sym)
  Matrix::forceSymmetric(dist, uplo = "U") %>%
    as.dist()
}

distance_frechet <- function(df, groups, by, ...) {
  groups_sym <- enquo(groups)
  by_sym <- enquo(by)
  vars_quo  <- enquos(...)
  l <- df %>%
    arrange(!!groups_sym, !!by_sym) %>%
    group_nest(!!groups_sym) %>%
    mutate(
      data = map(
        data,
        ~dplyr::select(.x, !!!vars_quo) %>%
          as.matrix()
      )
    )
  n <- nrow(l)
  dist <- matrix(NA_real_, n, n)
  plan(multisession(workers = 10))
  dist_df <- combn(seq_len(n), 2) %>%
    t() %>%
    magrittr::set_colnames(c("i", "j")) %>%
    as_tibble() %>%
    mutate(
      dist = future_map2_dbl(i, j, ~Frechet(l$data[[.x]], l$data[[.y]]), .progress = TRUE)
    )
  dist[dist_df$i, dist_df$j] <- dist_df$dist
  diag(dist) <- 0
  # for (i in seq_len(n)) {
  #   for (j in seq(i, n)) {
  #     dist[i, j] <- Frechet(l$data[[i]], l$data[[j]])
  #   }
  # }
  rownames(dist) <- pull(l, !!groups_sym)
  colnames(dist) <- pull(l, !!groups_sym)
  # browser()
  dist
  # Matrix::forceSymmetric(dist, uplo = "U") %>%
  #   as.dist()
}


# Calculate distances ----------------------------------------------------------
###############################################################################T

pc_lfc_long_z <- pc_lfc_long %>%
  standardize_data(gene_id, log2FoldChange, PC1, PC2) %>%
  inner_join(condition_meta, by = "condition")

pc_lfc_dist_euclidian <- distance_euclidian(
  pc_lfc_long_z %>%
    filter(Time == 24),
  gene_id, condition, PC1, PC2, log2FoldChange
)

# pc_lfc_dist_frechet <- distance_frechet(
#   pc_lfc_long_z %>%
#     filter(Time == 24),
#   gene_id, condition, PC1, PC2, log2FoldChange
# )
# write_rds(pc_lfc_dist_frechet, file.path(wd, "pc_lfc_frechet_dist.rds"))

# pc_lfc_dist_frechet <- read_rds(file.path(wd, "pc_lfc_frechet_dist.rds"))

# Do K-medoids -----------------------------------------------------------------
###############################################################################T

set.seed(42)
cluster_euclidian_raw <- pam(pc_lfc_dist_euclidian, k = 8, diss = TRUE)

cluster_euclidian_raw_df <- cluster_euclidian_raw %>%
  chuck("clustering") %>%
  enframe("gene_id", "cluster_id")

plot_cluster_trajectories(
  cluster_euclidian_raw_df %>%
    inner_join(
      pc_lfc_long_z %>%
        filter(Time == 24),
      by = "gene_id"
    ) %>%
    gather("PC", "PC_value", PC1, PC2) %>%
    inner_join(select(pca_values, condition, PC1), by = "condition"),
  PC_value, log2FoldChange, condition, PC1, gene_id, PC, cluster_id,
  all_traces = TRUE
)

cluster_euclidian_identities <- tribble(
  ~cluster_id, ~cluster, ~direction,
  1, "high_erk", "-",
  2, "low_erk", "-",
  3, "full_range", "+",
  4, "bell", "+",
  5, "low_erk", "+",
  6, "bell", "-",
  7, "high_erk", "+",
  8, "full_range", "-"
) %>%
  mutate(cluster_combined = paste(cluster, direction, sep = "_"))

cluster_euclidian <- cluster_euclidian_raw_df %>%
  inner_join(cluster_euclidian_identities, by = "cluster_id") %>%
  select(-cluster_id)

write_csv(
  cluster_euclidian,
  file.path(wd, "kmedoids_clustering_classes.csv")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T


activity <- Activity(
  "Clustering ERK response genes based on their expression profile with increasing ERK concentration",
  used = c(
    "syn21432183",
    "syn21432975",
    "syn21444456",
    "syn21444486"
  ),
  executed = "https://github.com/clemenshug/erk_senescence/blob/master/analysis/cluster_kmedoids.R"
)

syn_kmedoids_clusters <- Folder("kmedoids_clusters", parent = syn_clustering) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(wd, "kmedoids_clustering_classes.csv")
) %>%
  synStoreMany(syn_kmedoids_clusters, activity = activity)

