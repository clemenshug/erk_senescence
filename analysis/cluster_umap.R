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

wd <- here("umap")
dir.create(wd, recursive = TRUE, showWarnings = FALSE)

# set directories, import files ------------------------------------------------
###############################################################################T

meta <- syn("syn21432975") %>%
  read_csv()

cluster_overlap <- syn("syn21576614") %>%
  read_csv()

norm_counts <- syn("syn22093702") %>%
  read_csv()

cluster_names <- syn("syn21567677") %>%
  read_csv() %>%
  bind_rows(list(class_combined = "unknown", class_name = "Unknown")) %>%
  mutate_at(vars(class_name), . %>% as.factor() %>% fct_inorder())

temporal_ordering <- syn("syn21536903") %>%
  read_csv()

pca_loadings <- syn("syn21444455") %>%
  read_csv()

# UMAP embedding of variance-stabilized counts ---------------------------------
###############################################################################T

## All samples

vs_matrix <- norm_counts %>%
  filter(type == "normalized") %>%
  dplyr::select(ensembl_gene_id, sample_id, count) %>%
  spread(sample_id, count) %>%
  filter(ensembl_gene_id %in% function_clusters$gene_id) %>%
  column_to_rownames("ensembl_gene_id") %>%
  as.matrix()

vs_embedding_raw <- umap(
  vs_matrix, method = "naive",
  n_neighbors = 30, metric = "pearson"
)
vs_embedding <- vs_embedding_raw %>%
  chuck("layout") %>%
  as_tibble(rownames = "gene_id") %>%
  inner_join(
    cluster_overlap, by = "gene_id"
  ) %>%
  dplyr::select(-consensus_n, -gene_name) %>%
  gather("cluster_method", "cluster", -gene_id, -V1, -V2)

umap_plot <- ggplot(vs_embedding, aes(V1, V2, color = cluster)) +
  geom_point(alpha = 0.7, shape = 16) +
  facet_grid(vars(cluster_method, ))

ggsave(
  file.path(wd, "umap_all_samples.pdf"),
  umap_plot,
  width = 10, height = 8
)

vs_embedding_pca <- vs_embedding_raw %>%
  chuck("layout") %>%
  as_tibble(rownames = "gene_id") %>%
  inner_join(
    pca_loadings %>%
      filter(PC %in% c("PC1", "PC2")),
    by = "gene_id"
  ) %>%
  arrange(abs(value))

umap_plot <- ggplot(vs_embedding_pca, aes(V1, V2, color = value)) +
  geom_point(alpha = 0.7, shape = 16) +
  facet_wrap(vars(PC)) +
  scale_color_distiller(
    palette = "RdBu", limits = c(-0.06, 0.06), oob = scales::squish
  ) +
  labs(color = "PCA loading")

ggsave(
  file.path(wd, "umap_all_samples_pca.pdf"),
  umap_plot,
  width = 10, height = 6
)

# 24h samples

vs_matrix <- norm_counts %>%
  filter(
    type == "normalized",
    sample_id %in% {
      meta %>%
        filter(Time == 24) %>%
        pull(Sample_ID)
    }
  ) %>%
  dplyr::select(ensembl_gene_id, sample_id, count) %>%
  spread(sample_id, count) %>%
  filter(ensembl_gene_id %in% function_clusters$gene_id) %>%
  column_to_rownames("ensembl_gene_id") %>%
  as.matrix()

vs_embedding_raw <- umap(
  vs_matrix, method = "naive",
  n_neighbors = 30, metric = "pearson"
)
vs_embedding <- vs_embedding_raw %>%
  chuck("layout") %>%
  as_tibble(rownames = "gene_id") %>%
  inner_join(
    cluster_overlap, by = "gene_id"
  ) %>%
  dplyr::select(-consensus_n, -gene_name) %>%
  gather("cluster_method", "cluster", -gene_id, -V1, -V2)

umap_plot <- ggplot(vs_embedding, aes(V1, V2, color = cluster)) +
  geom_point(alpha = 0.7, shape = 16) +
  facet_wrap(vars(cluster_method))

ggsave(
  file.path(wd, "umap_24h_only.pdf"),
  width = 10, height = 8
)

vs_embedding_pca <- vs_embedding_raw %>%
  chuck("layout") %>%
  as_tibble(rownames = "gene_id") %>%
  inner_join(
    pca_loadings %>%
      filter(PC %in% c("PC1", "PC2")),
    by = "gene_id"
  ) %>%
  arrange(abs(value))

umap_plot <- ggplot(vs_embedding_pca, aes(V1, V2, color = value)) +
  geom_point(alpha = 0.7, shape = 16) +
  facet_wrap(vars(PC)) +
  scale_color_distiller(
    palette = "RdBu", limits = c(-0.06, 0.06), oob = scales::squish
  ) +
  labs(color = "PCA loading")

ggsave(
  file.path(wd, "umap_24h_only_pca.pdf"),
  umap_plot,
  width = 10, height = 6
)



