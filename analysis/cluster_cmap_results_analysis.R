library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(RColorBrewer)
library(pheatmap)
library(viridis)
library(ComplexHeatmap)
library(withr)
library(mice)
library(seriation)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("functional_enrichment")
dir.create(wd, recursive = TRUE, showWarnings = FALSE)

# set directories, import files ------------------------------------------------
###############################################################################T

function_clusters <- syn("syn21576614") %>%
  read_csv()

cluster_names <- syn("syn21567677") %>%
  read_csv() %>%
  bind_rows(list(class_combined = "unknown", class_name = "Unknown")) %>%
  mutate_at(vars(class_name), . %>% as.factor() %>% fct_inorder())

meta <- syn("syn21432975") %>%
  read_csv()

condition_meta <- meta %>%
  distinct(condition, DMSO, ERKi, Time, DOX)

cmap_results_all_conditions <- syn("syn22976182") %>%
  read_csv(col_types = "cccicdcccc")

cmap_results_consensus_clusters <- syn("syn22976162") %>%
  read_csv(col_types = "cccicdcccc")

pertubation_meta <- syn("syn21547097") %>%
  read_csv()

pca_values <- syn("syn21444456") %>%

# Plotting CMap heatmaps -------------------------------------------------------
###############################################################################T

cmap_results_cleaned <- cmap_results_consensus_clusters %>%
  mutate(
    col_id = paste(gene_set, cell_id, sep = "_")
  )

pertubation_meta_cleaned <- pertubation_meta %>%
  distinct(pert_id, pert_iname, pert_type) %>%
  bind_rows(
    cmap_results_consensus_clusters %>%
      filter(result_type == "pcl") %>%
      mutate(pert_id = id, pert_iname = id, pert_type = "pcl") %>%
      distinct(pert_id, pert_iname, pert_type)
  )

hcluster <- function(mat, dist_func = function(x) 1 - cor(x), linkage = "average") {
  dist_mat <- dist_func(mat) %>%
    as.dist()
  # browser()
  hclust(dist_mat, method = linkage) %>%
    seriation:::reorder.hclust(dist = dist_mat, method = "OLO")
}

plot_cmap_heatmap <- function(
  data,
  pert_filter_fun = . %>%
    mutate(direction = if_else(tau > 0, "correlated", "anti-correlated")) %>%
    arrange(desc(abs(tau))) %>%
    # group_by(gene_set, direction) %>%
    group_by(col_id) %>%
    slice(1:25) %>%
    pull(pert_id),
  row_order = NULL,
  col_order = NULL,
  ...
) {
  # browser()
  hm_data <- data %>%
    select(col_id, tau, pert_id) %>%
    # group_by(col_id, pert_id) %>%
    # summarize(tau = mean(tau), .groups = "drop") %>%
    filter(
      if (!is.null(pert_filter_fun))
        pert_id %in% pert_filter_fun(.)
      else TRUE
    ) %>%
    pivot_wider(id_cols = pert_id, names_from = col_id, values_from = tau)

  cmap_colormap <- circlize::colorRamp2(
    c(-100, -80, 80, 100),
    colors = c("#0000FF", "#FFFFFF", "#FFFFFF", "#FF0000"),
    space = "RGB"
  )

  row_annotation_data <- pertubation_meta_cleaned %>%
    mutate(
      name = paste(pert_iname, pert_type, sep = "_")
    )

  col_annotation_data <- data %>%
    distinct(col_id, gene_set, cell_id)

  hm_matrix <- hm_data %>%
    left_join(
      select(row_annotation_data, pert_id, name),
      by = "pert_id"
    ) %>%
    select(-pert_id) %>%
    # In rare cases two measurements for same pertubations, averaging
    group_by(name) %>%
    {
      d <- group_data(.)
      l <- map_int(d[[".rows"]], length)
      if (any(l > 1))
        warning("Pertubations with more than one measurement:", paste(d[["name"]][l > 1], collapse = ", "))
      summarize(., across(.fns = mean))
    } %>%
    column_to_rownames("name") %>%
    as.matrix()

  hm_matrix_imputed <- hm_matrix %>%
    mice() %>%
    mice::complete()

  # browser()

  if (is.null(col_order))
    col_cluster <- hcluster(hm_matrix_imputed)
  else {
    col_cluster <- FALSE
    hm_matrix <- hm_matrix[, col_order]
  }
  if (is.null(row_order))
    row_cluster <- hcluster(t(hm_matrix_imputed))
  else {
    row_cluster = FALSE
    hm_matrix <- hm_matrix[col_order, ]
  }

  # browser()

  row_annotation_data <- row_annotation_data %>%
    slice(match(rownames(hm_matrix), name))

  col_annotation_data <- col_annotation_data %>%
    slice(match(colnames(hm_matrix), col_id)) %>%
    inner_join(
      select(condition_meta, condition, Time, ERKi, DOX),
      by = c("gene_set" = "condition")
    ) %>%
    select(cell_id, Time, ERKi, DOX)

  Heatmap(
    hm_matrix,
    cluster_rows = row_cluster,
    cluster_columns = col_cluster,
    # row_dend_reorder = TRUE,
    # clustering_distance_rows = "pearson",
    right_annotation = rowAnnotation(
      `Pertubation type` = row_annotation_data[["pert_type"]]
    ),
    top_annotation = columnAnnotation(
      df = col_annotation_data %>%
        as.data.frame()
    ),
    col = cmap_colormap,
    name = "Tau",
    ...
  )
}

heatmap_pcl_cp_oe_sh_top_per_condition <- plot_cmap_heatmap(
  cmap_results_cleaned %>%
    filter(
      cell_id %in% c("A375", "summary"),
      pert_type %in% c("pcl", "trt_cp", "trt_oe", "trt_sh.cgs")
    ),
  . %>%
    mutate(direction = if_else(tau > 0, "correlated", "anti-correlated")) %>%
    arrange(desc(abs(tau))) %>%
    # group_by(gene_set, direction) %>%
    group_by(col_id) %>%
    slice(1:10) %>%
    pull(pert_id)
)

cmap_results_all_conditions_cleaned <- cmap_results_all_conditions %>%
  mutate(
    col_id = paste(gene_set, cell_id, sep = "_")
  )

heatmap_pcl_cp_oe_sh_top10_per_tp_24h_all <- cmap_results_all_conditions_cleaned %>%
  filter(
    gene_set %in% {
      condition_meta %>%
        filter(Time == 24) %>%
        pull(condition)
    },
    cell_id %in% c("A375", "summary"),
    pert_type %in% c("pcl", "trt_cp", "trt_oe", "trt_sh.cgs")
  ) %>%
  plot_cmap_heatmap(
    . %>%
      mutate(direction = if_else(tau > 0, "correlated", "anti-correlated")) %>%
      arrange(desc(abs(tau))) %>%
      # group_by(gene_set, direction) %>%
      group_by(col_id) %>%
      slice(1:10) %>%
      pull(pert_id),
    col_order = distinct(., condition = gene_set, cell_id, col_id) %>%
      inner_join(select(pca_values, condition, PC1), by = "condition") %>%
      arrange(cell_id, PC1) %>%
      pull(col_id)
  )

with_pdf(
  file.path(wd, "cmap_results_heatmap_pcl_cp_oe_sh_top10_per_tp_24h_all.pdf"),
  width = 8,
  height = 17,
  print(heatmap_pcl_cp_oe_sh_top10_per_tp_24h_all)
)


heatmap_pcl_cp_oe_sh_top10_dox1_time24 <- cmap_results_all_conditions_cleaned %>%
  filter(
    gene_set %in% {
      condition_meta %>%
        filter(Time == 24, DOX == 1) %>%
        pull(condition)
    },
    cell_id %in% c("A375", "summary"),
    pert_type %in% c("pcl", "trt_cp", "trt_oe", "trt_sh.cgs")
  ) %>%
  plot_cmap_heatmap(
    . %>%
      mutate(direction = if_else(tau > 0, "correlated", "anti-correlated")) %>%
      arrange(desc(abs(tau))) %>%
      # group_by(gene_set, direction) %>%
      group_by(col_id) %>%
      slice(1:10) %>%
      pull(pert_id),
    col_order = distinct(., condition = gene_set, cell_id, col_id) %>%
      inner_join(select(pca_values, condition, PC1), by = "condition") %>%
      arrange(cell_id, PC1) %>%
      pull(col_id)
  )

with_pdf(
  file.path(wd, "cmap_results_heatmap_pcl_cp_oe_sh_top10_dox1_time24.pdf"),
  width = 8,
  height = 17,
  print(heatmap_pcl_cp_oe_sh_top10_dox1_time24)
)

heatmap_pcl_cp_oe_sh_top10_dox0_time24 <- cmap_results_all_conditions_cleaned %>%
  filter(
    gene_set %in% {
      condition_meta %>%
        filter(Time == 24, DOX == 0) %>%
        pull(condition)
    },
    cell_id %in% c("A375", "summary"),
    pert_type %in% c("pcl", "trt_cp", "trt_oe", "trt_sh.cgs")
  ) %>%
  plot_cmap_heatmap(
    . %>%
      mutate(direction = if_else(tau > 0, "correlated", "anti-correlated")) %>%
      arrange(desc(abs(tau))) %>%
      # group_by(gene_set, direction) %>%
      group_by(col_id) %>%
      slice(1:10) %>%
      pull(pert_id),
    col_order = distinct(., condition = gene_set, cell_id, col_id) %>%
      inner_join(select(pca_values, condition, PC1), by = "condition") %>%
      arrange(cell_id, PC1) %>%
      pull(col_id)
  )

with_pdf(
  file.path(wd, "cmap_results_heatmap_pcl_cp_oe_sh_top10_dox0_time24.pdf"),
  width = 8,
  height = 17,
  print(heatmap_pcl_cp_oe_sh_top10_dox0_time24)
)


heatmap_pcl_cp_oe_sh_top10_dox1_ERK0 <- cmap_results_all_conditions_cleaned %>%
  filter(
    gene_set %in% {
      condition_meta %>%
        filter(DOX == 1, ERKi == 0) %>%
        pull(condition)
    },
    cell_id %in% c("A375", "summary"),
    pert_type %in% c("pcl", "trt_cp", "trt_oe", "trt_sh.cgs")
  ) %>%
  plot_cmap_heatmap(
    . %>%
      mutate(direction = if_else(tau > 0, "correlated", "anti-correlated")) %>%
      arrange(desc(abs(tau))) %>%
      # group_by(gene_set, direction) %>%
      group_by(col_id) %>%
      slice(1:10) %>%
      pull(pert_id),
    col_order = distinct(., condition = gene_set, cell_id, col_id) %>%
      inner_join(select(pca_values, condition, PC1), by = "condition") %>%
      inner_join(condition_meta, by = "condition") %>%
      arrange(cell_id, Time) %>%
      pull(col_id)
  )

with_pdf(
  file.path(wd, "cmap_results_heatmap_pcl_cp_oe_sh_top10_dox1_ERK0.pdf"),
  width = 8,
  height = 17,
  print(heatmap_pcl_cp_oe_sh_top10_dox1_ERK0)
)

heatmap_pcl_cp_oe_sh_top10_dox1_ERK1000 <- cmap_results_all_conditions_cleaned %>%
  filter(
    gene_set %in% {
      condition_meta %>%
        filter(DOX == 1, ERKi == 1000) %>%
        pull(condition)
    },
    cell_id %in% c("A375", "summary"),
    pert_type %in% c("pcl", "trt_cp", "trt_oe", "trt_sh.cgs")
  ) %>%
  plot_cmap_heatmap(
    . %>%
      mutate(direction = if_else(tau > 0, "correlated", "anti-correlated")) %>%
      arrange(desc(abs(tau))) %>%
      # group_by(gene_set, direction) %>%
      group_by(col_id) %>%
      slice(1:10) %>%
      pull(pert_id),
    col_order = distinct(., condition = gene_set, cell_id, col_id) %>%
      inner_join(select(pca_values, condition, PC1), by = "condition") %>%
      inner_join(condition_meta, by = "condition") %>%
      arrange(cell_id, Time) %>%
      pull(col_id)
  )

with_pdf(
  file.path(wd, "cmap_results_hheatmap_pcl_cp_oe_sh_top10_dox1_ERK1000.pdf"),
  width = 8,
  height = 17,
  print(heatmap_pcl_cp_oe_sh_top10_dox1_ERK1000)
)


heatmap_pcl_cp_oe_sh_top10_all_lol <- cmap_results_all_conditions_cleaned %>%
  filter(
    gene_set %in% {
      condition_meta %>%
        filter(DOX == 1) %>%
        pull(condition)
    },
    cell_id %in% c("A375", "summary"),
    pert_type %in% c("pcl", "trt_cp", "trt_oe", "trt_sh.cgs")
  ) %>%
  plot_cmap_heatmap(
    . %>%
      mutate(direction = if_else(tau > 0, "correlated", "anti-correlated")) %>%
      arrange(desc(abs(tau))) %>%
      # group_by(gene_set, direction) %>%
      group_by(col_id) %>%
      slice(1:10) %>%
      pull(pert_id),
    col_order = distinct(., condition = gene_set, cell_id, col_id) %>%
      # inner_join(select(pca_values, condition, PC1), by = "condition") %>%
      inner_join(condition_meta, by = "condition") %>%
      inner_join(
        pca_values %>%
          select(condition, PC1) %>%
          inner_join(condition_meta, by = "condition") %>%
          filter(Time == 24) %>%
          select(ERKi, DOX, PC1),
        by = c("ERKi", "DOX")
      ) %>%
      arrange(cell_id, PC1, Time) %>%
      pull(col_id),
    column_split = distinct(., condition = gene_set, cell_id, col_id) %>%
      # inner_join(select(pca_values, condition, PC1), by = "condition") %>%
      inner_join(condition_meta, by = "condition") %>%
      inner_join(
        pca_values %>%
          select(condition, PC1) %>%
          inner_join(condition_meta, by = "condition") %>%
          filter(Time == 24) %>%
          select(ERKi, DOX, PC1),
        by = c("ERKi", "DOX")
      ) %>%
      arrange(cell_id, PC1, Time) %>%
      mutate(x = paste(cell_id, PC1)) %>%
      pull(x),
    column_gap = unit(0, "mm"), border = TRUE
  )


with_pdf(
  file.path(wd, "cmap_results_hheatmap_pcl_cp_oe_sh_top10_all_lol.pdf"),
  width = 15,
  height = 50,
  print(heatmap_pcl_cp_oe_sh_top10_all_lol)
)
