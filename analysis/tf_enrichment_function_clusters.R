library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(enrichR)
library(RColorBrewer)
library(pheatmap)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("tf_enrichment", "function_clustering")
dir.create(wd, recursive = TRUE, showWarnings = FALSE)

# set directories, import files ------------------------------------------------
###############################################################################T

function_clusters <- syn("syn21478380") %>%
  read_csv()
meta <- syn("syn21432975") %>%
  read_csv()
surface_fit <- syn("syn21444486") %>%
  read_csv()

condition_meta <- meta %>%
  distinct(condition, ERKi, Time, DOX)

# Find TF enrichment for each class --------------------------------------------
###############################################################################T

selected_tf_dbs <- c(
  "TRRUST_Transcription_Factors_2019",
  "TRANSFAC_and_JASPAR_PWMs",
  "TF_Perturbations_Followed_by_Expression",
  "TF-LOF_Expression_from_GEO",
  "ENCODE_TF_ChIP-seq_2015",
  "ChEA_2016"
  # "ARCHS4_TFs_Coexp"
)

enrichr_res_raw <- function_clusters %>%
  group_nest(class, direction, class_combined) %>%
  mutate(
    data = map(
      data,
      ~enrichr(.x$gene_name, databases = selected_tf_dbs)
    )
  )

enrichr_res <- enrichr_res_raw %>%
  mutate(
    data = map(
      data,
      enframe, "db", "data"
    )
  ) %>%
  unnest(data) %>%
  mutate_at(vars(data), map, . %>% as_tibble %>% arrange(desc(Combined.Score)))

# Plot TF enrichment as heatmap ------------------------------------------------
###############################################################################T


plot_heatmap <- function(mat, ...) {
  abs_max <- as.vector(mat) %>%
    abs() %>%
    quantile(.90, names = FALSE, na.rm = TRUE) %>%
    round(1)
  breaks <- seq(-abs_max, abs_max, by = 0.1)
  cmap <- rev(colorRampPalette(brewer.pal(7, "RdBu"))(length(breaks)))
  # row_clust <- hclust(as.dist(1 - cor(t(mat), method =  "pearson")), "average")
  row_clust <- hclust(dist(mat, method = "euclidian"), "ward.D2")
  # browser()
  pheatmap(
    mat,
    color = cmap,
    breaks = breaks,
    cluster_rows = row_clust,
    cluster_cols = FALSE,
    silent = TRUE,
    ...
    # cutree_rows = 6
  )$gtable
}

min_p_val <- 0.01

enrichr_mats <- enrichr_res %>%
  unnest(data) %>%
  group_nest(db) %>%
  mutate(
    top_10_tfs = map(
      data,
      ~.x %>%
        group_by(class_combined) %>%
        arrange(Adjusted.P.value, .by_group = TRUE) %>%
        mutate(rank = 1:n()) %>%
        ungroup() %>%
        filter(Adjusted.P.value <= min_p_val, rank <= 10) %>%
        distinct(Term, Adjusted.P.value)  %>%
        arrange(Adjusted.P.value) %>%
        pull(Term)
    ),
    passing_tfs = map(
      data,
      ~.x %>%
        group_by(Term) %>%
        arrange(Adjusted.P.value, .by_group = TRUE) %>%
        mutate(rank = 1:n()) %>%
        ungroup() %>%
        filter(Adjusted.P.value <= min_p_val, rank == 1) %>%
        distinct(Term, Adjusted.P.value) %>%
        arrange(Adjusted.P.value) %>%
        pull(Term)
    ),
    selected_tfs = map2(
      top_10_tfs, passing_tfs,
      ~c(.x, .y) %>%
        head(n = 50)
    ),
    mat = map2(
      data, selected_tfs,
      ~filter(.x, Term %in% .y) %>%
        select(Term, Combined.Score, class_combined) %>%
        spread(class_combined, Combined.Score, fill = 0) %>%
        column_to_rownames("Term") %>%
        as.matrix()
    )
  )

enrichr_heatmaps <- enrichr_mats %>%
  mutate(
    plot = map2(
      mat, data,
      ~plot_heatmap(
        .x,
        annotation_col = distinct(.y, class_combined, class, direction) %>%
          column_to_rownames("class_combined"),
        annotation_colors = list(
          direction = c("-" = "blue", "+" = "red", "0" = "grey"),
          class = brewer.pal(5, "Set2") %>%
            set_names(unique(.y$class))
        )
      )
    )
  )

pwalk(
  enrichr_heatmaps,
  function(db, plot, ...) {
    cowplot::ggsave2(
      file.path(wd, paste0("tf_enrichment_heatmap_top_10_per_class_", db, ".pdf")),
      plot
    )
  }
)
