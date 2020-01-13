library(tidyverse)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

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


# Temporal clustering ----------------------------------------------------------
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
        mat <- matrix(NA_real_, n, n, dimnames = list(d$gene_name, d$gene_name))
        for (i in 1:n) {
          for (j in i:n) {
            mat[[j, i]] <- sqrt(sum((d[["data"]][[i]][["log2FoldChange"]] - d[["data"]][[j]][["log2FoldChange"]])**2))
          }
        }
        as.dist(mat, diag = TRUE)
      }
    )
  )

temporal_clusters <- temporal_distance_per_cluster %>%
  mutate(
    data = map(
      data,
      hclust,
      method = "average"
    )
  )

plot_single_surface <- function(d) {
  abs_max <- quantile(abs(d$log2FoldChange), .95)
  d %>%
    mutate_at(vars(Time, ERKi), as.factor) %>%
    ggplot(aes(ERKi, Time, fill = log2FoldChange)) +
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

plot_single_line_by_erki(
  lfc_long %>%
    inner_join(condition_meta, by = "condition") %>%
    filter(gene_name == "JUNB")
)

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

lfc_long_z %>%
  inner_join(function_clusters, by = c("gene_id", "gene_name")) %>%
  plot_clusters_line_by_erki(class_combined)

# Temporal clustering based only on most extreme condition ---------------------
###############################################################################T



# First find condition with most extreme change at 24h

extreme_condition <- lfc_long %>%
  left_join(condition_meta, by = "condition") %>%
  filter(Time == 24, DOX == 1, gene_id %in% function_clusters$gene_id) %>%
  mutate(log2FoldChange_abs = abs(log2FoldChange)) %>%
  arrange(gene_id, gene_name, desc(log2FoldChange_abs)) %>%
  group_by(gene_id, gene_name) %>%
  slice(1) %>%
  ungroup()

# Then find time trajectory for each gene for that condition

time_course <- lfc_long %>%
  left_join(condition_meta, by = "condition") %>%
  semi_join(extreme_condition, by = c("gene_id", "ERKi", "DOX")) %>%
  arrange(gene_name, gene_id, Time)

# Add time induction metrics, like 50% induction, extremum etc

induction_metrics <- time_course %>%
  group_by(gene_id, gene_name) %>%
  summarize(
    half_induction =
  )

time_course %>%
  filter(gene_id %in% sample(unique(.$gene_id), 20)) %>%
  ggplot(aes(Time, log2FoldChange)) +
  facet_wrap(vars(gene_name)) +
  geom_line()

time_course_data <- lfc_long %>%
  inner_join(
    filter(time_course_meta, ERKi %in% c(0, 1000)),
    # time_course_meta,
    by = "condition"
  ) %>%
  group_by(gene_id) %>%
  mutate(log2FoldChange_scaled = scale(log2FoldChange), log2FoldChange_sdev = sd(log2FoldChange)) %>%
  ungroup() %>%
  inner_join(select(pc_clusters, gene_id, cluster), by = "gene_id") %>%
  group_by(ERKi, cluster, gene_id) %>%
  ungroup() %>%
  group_nest(ERKi, cluster)
  # group_by(ERKi) %>%
  # mutate(log2FoldChange_sdev = bind_rows(data) %>% pull(log2FoldChange) %>% sd()) %>%
  # ungroup()

time_course_matrix <- time_course_data %>%
  mutate(
    log2FoldChange_mat = map(
      data,
      ~select(.x, gene_id, Time, log2FoldChange_scaled) %>%
        spread(Time, log2FoldChange_scaled) %>%
        column_to_rownames("gene_id") %>%
        as.matrix()
    ),
    sdev_cross_matrix = map(
      data,
      ~group_by(.x, gene_id) %>%
        summarize(log2FoldChange_sdev = sd(log2FoldChange)) %>%
        ungroup() %>%
        {set_names(.$log2FoldChange_sdev, .$gene_id)} %>%
        t() %>%
        crossprod()
    )
  )

# time_course_dist_matrix <- time_course_matrix %>%
#   mutate(
#     data = map(
#       data,
#       ~t(1 - cor(t(.x), method = "spearman"))
#     )
#   )
time_course_dist_matrix <- time_course_matrix %>%
  mutate(
    data = map(
      log2FoldChange_mat,
      ~dist(.x, method = "euclidian")
    )
  )

time_course_dist_matrix_agg <- time_course_dist_matrix %>%
  group_by(cluster) %>%
  summarize(
    data = reduce(data, magrittr::add) %>%
      list()
    # data = map(data, as.matrix) %>%
    #   map2(sdev_cross_matrix, magrittr::multiply_by) %>%
    #   map(as.dist) %>%
    #   reduce(magrittr::add) %>%
    #   # magrittr::multiply_by(1 / sum(log2FoldChange_sdev)) %>%
    #   list()
  ) %>%
  ungroup()

time_course_clust <- time_course_dist_matrix %>%
  mutate(
    data = map(
      data,
      ~cluster::pam(.x, k = 5)
    )
  )

time_course_clust_df <-  time_course_clust %>%
  mutate(
    data = map(
      data,
      ~chuck(.x, "clustering") %>%
        enframe("gene_id", "cluster") %>%
        mutate_at(vars(cluster), as.character)
    )
  )

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

time_course_clust_traces <- time_course_clust_df %>%
  left_join(rename(time_course_data, tc_data = data), by = c("cluster", "ERKi")) %>%
  mutate(
    data = map2(
      data, tc_data,
      left_join, by = "gene_id"
    ) %>%
      map(rename, tc_cluster = cluster)
  ) %>%
  rename(pc_cluster = cluster) %>%
  select(-tc_data) %>%
  unnest(data)


plot_cluster_trajectories(
  time_course_clust_traces %>%
    filter(pc_cluster == 2) %>%
    mutate(log2FoldChange_normed = log2FoldChange_scaled),
  Time, log2FoldChange_normed, Time, gene_id, ERKi, tc_cluster,
  all_traces = TRUE
)

time_course_clust_hm <- time_course_clust_df %>%
  mutate(
    plot = map(
      data,
      function(df) {
        mat <- lfc_long %>%
          filter(gene_id %in% df$gene_id) %>%
          select(gene_id, )
      }
    )
  )
