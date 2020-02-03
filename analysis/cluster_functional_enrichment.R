library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(enrichR)
library(RColorBrewer)
library(pheatmap)
library(ggforce)
library(viridis)
library(furrr)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("functional_enrichment")
dir.create(wd, recursive = TRUE, showWarnings = FALSE)

# set directories, import files ------------------------------------------------
###############################################################################T

function_clusters <- syn("syn21478380") %>%
  read_csv()

deseq_padj <- syn("syn21432184") %>%
  read_csv()

cluster_names <- syn("syn21567677") %>%
  read_csv() %>%
  mutate_at(vars(class_name), . %>% as.factor() %>% fct_inorder())

temporal_ordering <- syn("syn21536903") %>%
  read_csv()

# EnrichR ----------------------------------------------------------------------
###############################################################################T

function_clusters_enrichr <- function_clusters %>%
  group_by(class, direction, class_combined) %>%
  summarize(
    enrichment = list(
      enrichr(
        gene_name,
        databases = c(
          "GO_Biological_Process_2018",
          "GO_Molecular_Function_2018",
          "KEGG_2019_Human",
          "Reactome_2016"
        )
      ) %>%
        bind_rows(.id = "database") %>%
        as_tibble() %>%
        arrange(desc(Combined.Score))
    )
  ) %>%
  ungroup()

# TopGO functions --------------------------------------------------------------
###############################################################################T

library(topGO)
all_genes <- deseq_padj %>%
  pull(gene_id) %>%
  unique()

go_objects <- list(
  "bp" = new(
    getClassDef("topGOdata", package = "topGO"),
    ontology = "BP",
    allGenes = set_names(rep(1, length(all_genes)), all_genes),
    geneSelectionFun = function (x) x > 0.5,
    annot = topGO::annFUN.org,
    mapping = "org.Hs.eg.db",
    ID = "ensembl"
  ),
  "mf" = new(
    getClassDef("topGOdata", package = "topGO"),
    ontology = "MF",
    allGenes = set_names(rep(1, length(all_genes)), all_genes),
    geneSelectionFun = function (x) x > 0.5,
    annot = topGO::annFUN.org,
    mapping = "org.Hs.eg.db",
    ID = "ensembl"
  )
)




topgo_enrichment <- function(gene_set, all_genes, go_domain = "bp", ...) {
  gene_input <- set_names(
    if_else(all_genes %in% gene_set, 1, 0),
    all_genes
  )
  GOdata <- topGO::updateGenes(
    go_objects[[go_domain]], gene_input, function (x) x > 0.5
  )
  suppressMessages(resultFisher <- topGO::runTest(
    GOdata,
    algorithm = "weight01",
    statistic = "fisher"
  ))
  # resultKS <- topGO::runTest(
  #   GOdata,
  #   algorithm = "weight01",
  #   statistic = "ks"
  # )
  topGO::GenTable(
    GOdata,
    fisher = resultFisher,
    # ks = resultKS,
    orderBy = "fisher",
    # topNodes = max(length(resultFisher@score), length(resultKS@score)),
    topNodes = length(resultFisher@score),
    numChar = 1000
  ) %>%
    as_tibble() %>%
    dplyr::rename(id = GO.ID, term = Term, pval = fisher, annotated = Annotated, significant = Significant) %>%
    dplyr::mutate_at(vars(pval), ~as.numeric(gsub("< 1e-30", "1e-30", .x)))
}

# topGO on function clusters ---------------------------------------------------
###############################################################################T

function_clusters_topgo <- function_clusters %>%
  group_nest(class, direction, class_combined) %>%
  crossing(go_domain = c("bp", "mf")) %>%
  mutate(
    enrichment = map2(
      data, go_domain,
      ~topgo_enrichment(.x$gene_id, all_genes, go_domain = .y)
    )
  )


function_clusters_topgo_df <- function_clusters_topgo %>%
  dplyr::select(-data) %>%
  unnest(enrichment) %>%
  arrange(pval) %>%
  filter(
    id %in% (
        c(
        filter(., pval <= 0.05) %>%
          group_by(class, direction, class_combined) %>%
          dplyr::slice(1:5) %>%
          ungroup() %>%
          pull(id),
        filter(., pval <= 0.05) %>%
          pull(id)
      ) %>%
        unique() %>%
        head(50)
    )
  )

plot_heatmap <- function(mat, ...) {
  abs_max <- c(mat) %>%
    abs() %>%
    quantile(.95, names = FALSE, na.rm = TRUE) %>%
    round(1)
  breaks <- seq(0, abs_max, by = 0.1)
  cmap <- colorRampPalette(c("#ffffff", brewer.pal(7, "Reds")))(length(breaks))
  # cmap <- magma(length(breaks) - 1)
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

function_clusters_topgo_mat <- function_clusters_topgo_df %>%
  inner_join(cluster_names, by = "class_combined") %>%
  transmute(
    class_name,
    term = paste0(go_domain, "_", term),
    pval = -log10(pval)
  ) %>%
  spread(class_name, pval) %>%
  column_to_rownames("term") %>%
  as.matrix()

function_clusters_topgo_hm <- plot_heatmap(
  function_clusters_topgo_mat
    # magrittr::inset(. < -log10(0.05), 0)
)
ggsave(
  file.path(wd, "function_clusters_go_heatmap.pdf"),
  function_clusters_topgo_hm, width = 9, height = 12
)
grid::grid.draw(function_clusters_topgo_hm)

# topGO on time series induction -----------------------------------------------
###############################################################################T

temporal_ordering_abs <- temporal_ordering %>%
  filter(directed == "absolute") %>%
  mutate_at(
    vars(max_induction, mid_induction),
    cut, breaks = c(0, 3, 10, Inf), labels = c("early", "mid", "late")
  ) %>%
  mutate(
    max_induction = fct_cross(max_induction, direction_max, sep = "_"),
    mid_induction = fct_cross(mid_induction, direction_max, sep = "_"),
  )

plan(multisession(workers = 8))
temporal_ordering_abs_mid_ind_topgo <- temporal_ordering_abs %>%
  filter(gene_id %in% function_clusters$gene_id) %>%
  group_nest(mid_induction, ERKi) %>%
  crossing(go_domain = c("bp", "mf")) %>%
  mutate(
    enrichment = future_map2(
      data, go_domain,
      ~topgo_enrichment(.x$gene_id, all_genes, go_domain = .y),
      .progress = TRUE
    )
  )

