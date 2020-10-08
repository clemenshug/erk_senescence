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

# Plotting CMap heatmaps -------------------------------------------------------
###############################################################################T

cmap_results_cleaned <- cmap_results_consensus_clusters %>%
  mutate(
    pert_id = if_else(is.na(pert_id), id, pert_id),
    pert_type = if_else(is.na(pert_type), "pert_pcl", pert_type),
    cell_id = if_else(is.na(cell_id), "summary", cell_id),
    col_id = paste(gene_set, cell_id, sep = "_")
  ) %>%
  drop_na(pert_id)

pertubation_meta_cleaned <- pertubation_meta %>%
  distinct(pert_id, pert_iname, pert_type) %>%
  bind_rows(
    cmap_results_consensus_clusters %>%
      filter(result_type == "pcl") %>%
      mutate(pert_id = id, pert_iname = id, pert_type = "pert_pcl") %>%
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
    pull(pert_id)
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

  row_annotation_data <- pertubation_meta_cleaned %>%
    slice(match(hm_data[["pert_id"]], pert_id)) %>%
    mutate(
      name = paste(pert_iname, pert_type, sep = "_")
    )

  col_annotation_data <- data %>%
    distinct(col_id, gene_set, cell_id) %>%
    slice(match(colnames(hm_data), col_id))

  cmap_colormap <- circlize::colorRamp2(
    c(-100, -80, 80, 100),
    colors = c("#0000FF", "#FFFFFF", "#FFFFFF", "#FF0000"),
    space = "RGB"
  )
  # browser()

  hm_matrix <- hm_data %>%
    left_join(
      select(row_annotation_data, pert_id, name),
      by = "pert_id"
    ) %>%
    select(-pert_id) %>%
    column_to_rownames("name") %>%
    as.matrix()

  hm_matrix_imputed <- hm_matrix %>%
    mice() %>%
    mice::complete()

  # browser()

  col_cluster <- hcluster(hm_matrix_imputed)
  row_cluster <- hcluster(t(hm_matrix_imputed))

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
      `Cell type` = col_annotation_data[["cell_id"]]
    ),
    col = cmap_colormap,
    name = "Tau"
  )
}

heatmap_pcl_cp_oe_sh_top_per_condition <- plot_cmap_heatmap(
  cmap_results_cleaned %>%
    filter(
      cell_id %in% c("A375", "summary"),
      pert_type %in% c("pert_pcl", "trt_cp", "trt_oe", "trt_sh.cgs")
    )
)

with_pdf(
  file.path(wd, "cmap_results_heatmap__pcl_cp_oe_sh_top_per_condition.pdf"),
  width = 5,
  height = 15,
  print(hm)
)
