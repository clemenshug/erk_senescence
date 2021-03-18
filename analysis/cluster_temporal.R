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
surface_fit <- syn("syn21444486") %>%
  read_csv()

condition_meta <- meta %>%
  distinct(condition, ERKi, Time, DOX) %>%
  # Expand conditions to include explicitely the initial condition
  bind_rows(
    crossing(
      Time = 0,
      ERKi = unique(.[["ERKi"]]),
      DOX = c(0, 1)
    ) %>%
      mutate(condition = paste0("time", Time, "dox", DOX, "conc", ERKi))
  )

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
  # Add explicit zeroes for initial condition
  bind_rows(
    distinct(., gene_id, gene_name) %>%
      crossing(
        log2FoldChange = 0,
        condition = setdiff(condition_meta$condition, names(pairwise_lfc))
      )
  )

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
temporal_lfc <- lfc_long %>%
  inner_join(condition_meta, by = "condition") %>%
  filter(DOX == 1, gene_id %in% surface_fit$gene_id) %>%
  arrange(gene_id, gene_name, ERKi, Time) %>%
  select(-DOX)

temporal_var_erki <- temporal_lfc %>%
  group_by(gene_id, gene_name, ERKi) %>%
  summarize(
    log2FoldChange_var_erki = var(log2FoldChange),
    log2FoldChange_range_erki = max(log2FoldChange) - min(log2FoldChange),
    .groups = "drop"
  )

# Summarize LFC across all ERKi concentrations
temporal_lfc_ERKi_sum <- temporal_lfc %>%
  inner_join(
    temporal_var_erki,
    by = c("gene_id", "gene_name", "ERKi")
  ) %>%
  group_by(
    gene_id, gene_name, Time
  ) %>%
  summarize(
    # Weigh ERK concentrations by the gene's variance across time
    log2FoldChange_sum_mean = weighted.mean(log2FoldChange, log2FoldChange_range_erki),
    log2FoldChange_sum_var_erki = Hmisc::wtd.var(log2FoldChange, log2FoldChange_range_erki, method = "ML"),
    # Use either 33rd or 66th percentile, whichever absolute value is larger
    # CMap method for aggregating cell lines
    log2FoldChange_sum_cmap = Hmisc::wtd.quantile(
      log2FoldChange, probs = c(0.67, 0.33), normwt = TRUE,
      weights = log2FoldChange_range_erki
    ) %>%
      {.[order(abs(.))[2]]},
    ERKi = "all",
    .groups = "drop"
  ) %>%
  group_by(gene_id, gene_name) %>%
  mutate(
    log2FoldChange_sum_range_erki = max(log2FoldChange_sum_cmap) - min(log2FoldChange_sum_cmap)
  ) %>%
  ungroup() %>%
  select(
    gene_id, gene_name, Time, ERKi, log2FoldChange = log2FoldChange_sum_cmap,
    log2FoldChange_var_erki = log2FoldChange_sum_var_erki,
    log2FoldChange_range_erki = log2FoldChange_sum_range_erki
  )

temporal_lfc_combined <- bind_rows(
  temporal_lfc %>%
    inner_join(
      temporal_var_erki,
      by = c("gene_id", "gene_name", "ERKi")
    ) %>%
    mutate(ERKi = as.character(ERKi)) %>%
    select(-condition),
  temporal_lfc_ERKi_sum
)

# Scale LFC across time for each gene
temporal_lfc_rescaled <- temporal_lfc_combined %>%
  arrange(gene_id, gene_name, ERKi, Time) %>%
  group_by(gene_id, gene_name, ERKi) %>%
  transmute(
    direction_max = if_else(
      log2FoldChange[order(abs(log2FoldChange), decreasing = TRUE)[[1]]] > 0,
      "pos",
      "neg"
    ),
    direction_end = if_else(
      tail(log2FoldChange, 1) > 0,
      "pos",
      "neg"
    ),
    Time,
    log2FoldChange,
    directed = log2FoldChange %>%
      # magrittr::divide_by(log2FoldChange[order(abs(log2FoldChange), decreasing = TRUE)[1]]),
      magrittr::divide_by(max(abs(.))),
    absolute = log2FoldChange %>%
      abs() %>%
      magrittr::divide_by(max(.))
  ) %>%
  ungroup() %>%
  gather(key = "directed", value = "log2FoldChange_rescaled", directed, absolute)

# cor_across_erk <- temporal_stats %>%
#   group_by(gene_id, gene_name, directed) %>%
#   group_modify(
#     function(df, g) {
#       # browser()
#       lfc_mat <- df %>%
#         select(ERKi, Time, log2FoldChange_rescaled) %>%
#         spread(ERKi, log2FoldChange_rescaled) %>%
#         column_to_rownames("Time") %>%
#         as.matrix()
#       var_vec <- df %>%
#         distinct(ERKi, log2FoldChange_var_erki) %>%
#         column_to_rownames("ERKi") %>%
#         as.matrix()
#       lfc_cor_mat <- cor(lfc_mat, method = "pearson")
#       var_cross_mat <- var_vec %*% t(var_vec)
#       weighted_mean_cor <- weighted.mean(c(lfc_cor_mat), c(var_cross_mat))
#       # browser()
#       tibble(weighted_mean_cor = weighted_mean_cor)
#     }
#   ) %>%
#   ungroup()

temporal_ordering <- temporal_lfc_rescaled %>%
  group_by(gene_id, gene_name, ERKi, directed) %>%
  summarize(
    max_induction_idx = order(abs(log2FoldChange_rescaled), decreasing = TRUE)[1],
    max_induction_time = Time[max_induction_idx],
    max_induction_rescaled = log2FoldChange_rescaled[max_induction_idx],
    max_induction = log2FoldChange[max_induction_idx],
    mid_induction_time = {
      if (directed == "absolute" || max_induction > 0)
        Time[min(which(log2FoldChange_rescaled > 0.5))]
      else
        Time[min(which(log2FoldChange_rescaled < -0.5))]
    },
    mean_induction_rescaled = mean(log2FoldChange_rescaled),
    mean_induction = mean(log2FoldChange),
    .groups = "drop"
  )

# plot_single_surface(
#   lfc_long %>%
#     inner_join(condition_meta, by = "condition") %>%
#     filter(gene_name == "MYC")
# )
#
# plot_single_surface(
#   lfc_rescaled %>%
#     filter(gene_name == "MYC") %>%
#     mutate(log2FoldChange_rescaled_norm = log2FoldChange_rescaled*log2FoldChange_var_time*if_else(direction_max == "pos", 1, -1)),
#   aes(ERKi, Time, fill = log2FoldChange_rescaled_norm)
# )

hm_mid_vs_max <- temporal_ordering %>%
  count(directed, mid_induction_time, max_induction_time, ERKi) %>%
  mutate_at(vars(mid_induction_time, max_induction_time), . %>% as.character() %>%  fct_inseq()) %>%
  ggplot(aes(mid_induction_time, max_induction_time, fill = n)) +
    geom_raster() +
    facet_grid(vars(ERKi), vars(directed))
cowplot::ggsave2(
  file.path(wd, "temporal_ordering_heatmap_mid_vs_max.pdf"),
  hm_mid_vs_max, width = 5, height = 8
)

barplot_mid_vs_max <- temporal_ordering %>%
  mutate_at(vars(mid_induction_time, max_induction_time), . %>% as.character() %>%  fct_inseq()) %>%
  gather("key", "value", max_induction_time, mid_induction_time) %>%
  mutate_at(vars(value), . %>% as.character() %>%  fct_inseq()) %>%
  ggplot(aes(value, fill = key)) +
    geom_bar(position = "dodge") +
    facet_grid(vars(ERKi), vars(directed))
cowplot::ggsave2(
  file.path(wd, "temporal_ordering_barplot_mid_vs_max.pdf"),
  barplot_mid_vs_max, width = 4, height = 8
)

write_csv(
  temporal_lfc_combined,
  file.path(wd, "temporal_lfc.csv")
)

write_csv(
  temporal_lfc_rescaled,
  file.path(wd, "temporal_lfc_rescaled.csv")
)

write_csv(
  temporal_ordering,
  file.path(wd, "temporal_ordering.csv")
)

# Write special JY combined table with log2FoldChange --------------------------
###############################################################################T

# temporal_ordering_with_lfc <- temporal_ordering %>%
#   filter(ERKi != "all") %>%
#   mutate_at(vars(ERKi), as.numeric) %>%
#   inner_join(
#     function_clusters,
#     by = c("gene_id", "gene_name")
#   ) %>%
#   nest_join(
#     temporal_lfc %>%
#       transmute(
#         gene_id,
#         ERKi,
#         Time = factor(as.character(Time), levels = as.character(sort(unique(Time)))) %>%
#           fct_relabel(~paste0("lfc_", .x)),
#         log2FoldChange
#       ) %>%
#       arrange(Time) %>%
#       spread(Time, log2FoldChange),
#     by = c("gene_id", "ERKi"),
#     name = "lfc_data"
#   ) %>%
#   unnest(lfc_data) %>%
#   arrange(gene_name, directed, ERKi)
#
# write_csv(
#   temporal_ordering_with_lfc,
#   file.path(wd, "temporal_ordering_with_lfc.csv")
# )

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
  file.path(wd, "temporal_lfc.csv"),
  file.path(wd, "temporal_ordering.csv"),
  file.path(wd, "temporal_lfc_rescaled.csv")
  # file.path(wd_clustering, "temporal_euclidian_hclust.rds"),
  # file.path(wd_clustering, "temporal_kmedoids.rds")
) %>%
  synStoreMany(temporal_ordering_syn, activity = activity)
