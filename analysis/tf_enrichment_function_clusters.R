library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(enrichR)
library(enrichr.db)
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

# Aggregate ChEA TF data -------------------------------------------------------
###############################################################################T

chea_tfs <- enrichr.db::enrichr_terms %>%
  filter(library == "ChEA_2016") %>%
  chuck("data", 1) %>%
  enframe("tf_info", "targets") %>%
  mutate(
    tf_info_split = str_split(tf_info, fixed("_")),
    tf = map_chr(tf_info_split, 1),
    species = map_chr(tf_info_split, tail, n = 1) %>%
      str_to_lower()
  ) %>%
  filter(species %in% c("human", "mouse", "rat"))

chea_tfs_agg <- chea_tfs %>%
  group_by(tf) %>%
  filter(if ("human" %in% species) species == "human" else TRUE) %>%
  # summarize(targets = list(reduce(targets, union, .init = c()))) %>%
  ungroup()

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

# Aggregate ChEA result
###############################################################################T



enrichr_res_agg <- enrichr_res %>%
  bind_rows(
    enrichr_res %>%
      filter(db == "ChEA_2016") %>%
      mutate(
        db = "ChEA_2016_agg",
        data = map(
          data,
          function(data) {
            data %>%
              inner_join(chea_tfs_agg, by = c("Term" = "tf_info")) %>%
              group_by(tf) %>%
              summarize(
                Adjusted.P.value = psych::harmonic.mean(Adjusted.P.value),
                Combined.Score = median(Combined.Score)
              ) %>%
              ungroup() %>%
              mutate(Term = tf)
          }
        )
      )
  )

write_rds(
  enrichr_res_agg,
  file.path(wd, "enrichr_results_functional_clusters.rds")
)

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

enrichr_mats <- enrichr_res_agg %>%
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
      plot,
      width = 10, height = 14
    )
  }
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

activity <- Activity(
  "Enrichment of TF targets in functional clusters",
  used = c(
    "syn21478380",
    "syn21432975",
    "syn21444486"
  ),
  executed = "https://github.com/clemenshug/erk_senescence/blob/master/analysis/tf_enrichment_function_clusters.R"
)

syn_enrichr <- "syn21492378"

c(
  file.path(wd, "enrichr_results_functional_clusters.rds")
) %>%
  synStoreMany(syn_enrichr, activity = activity)

