library(tidyverse)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

# set directories, import files ------------------------------------------------
###############################################################################T

pc_clusters <- syn("syn21432143") %>%
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
  gather("condition", "log2FoldChange", -gene_name, -gene_id) %>%
  group_by(gene_id, gene_name, condition) %>%
  summarize(log2FoldChange = mean(log2FoldChange)) %>%
  ungroup()

# Temporal clustering ----------------------------------------------------------
###############################################################################T

time_course_meta <- condition_meta %>%
  filter(DOX == 1)

time_course_data <- lfc_long %>%
  inner_join(
    # filter(time_course_meta, ERKi %in% c(0, 1000)),
    time_course_meta,
    by = "condition"
  ) %>%
  inner_join(select(pc_clusters, gene_id, cluster), by = "gene_id") %>%
  group_by(ERKi) %>%
  group_by(ERKi, cluster, gene_id) %>%
  mutate(log2FoldChange_scaled = scale(log2FoldChange), log2FoldChange_sdev = sd(log2FoldChange)) %>%
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
    data = map(data, as.matrix) %>%
      map2(sdev_cross_matrix, magrittr::multiply_by) %>%
      map(as.dist) %>%
      reduce(magrittr::add) %>%
      # magrittr::multiply_by(1 / sum(log2FoldChange_sdev)) %>%
      list()
  ) %>%
  ungroup()

time_course_clust <- time_course_dist_matrix_agg %>%
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
  left_join(rename(time_course_data, tc_data = data), by = "cluster") %>%
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
    mutate(log2FoldChange_normed = log2FoldChange),
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
