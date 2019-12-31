library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(enrichR)
library(RColorBrewer)
library(pheatmap)
library(ggforce)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("clustering", "overlap")
dir.create(wd, recursive = TRUE, showWarnings = FALSE)

# set directories, import files ------------------------------------------------
###############################################################################T

km_clusters <- syn("syn21432143") %>%
  read_csv() %>%
  transmute(
    gene_id,
    gene_name,
    cluster = recode(grouped_cluster, high_range = "high_erk", low_range = "low_erk", linear = "full_range"),
    direction = recode(response, negative = "-", positive = "+"),
    cluster_combined = paste0(cluster, "_", direction)
  )
lm_clusters <- syn("syn21448572") %>%
  read_csv() %>%
  rename(cluster = class, cluster_combined = class_combined) %>%
  select(-low, -high)
function_clusters <- syn("syn21478380") %>%
  read_csv() %>%
  rename(cluster = class, cluster_combined = class_combined)
meta <- syn("syn21432975") %>%
  read_csv()
surface_fit <- syn("syn21444486") %>%
  read_csv()

# Merge clustering data --------------------------------------------------------
###############################################################################T

cluster_algos <- tribble(
  ~algorithm, ~clusters,
  "k_medoids", km_clusters,
  "linear_model", lm_clusters,
  "function_fitting", function_clusters
)

cluster_overlap <- cluster_algos %>%
  unnest(clusters) %>%
  select(algorithm, gene_id, gene_name, cluster_combined) %>%
  spread(algorithm, cluster_combined) %>%
  drop_na()

cluster_overlap_venn <- cluster_overlap %>%
  group_by_at(vars(3:5)) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  gather_set_data(1:3) %>%
  ggplot(aes(x, id = id, split = y, value = n)) +
    geom_parallel_sets(aes(fill = k_medoids), alpha = 0.3, axis.width = 0.1) +
    geom_parallel_sets_axes(axis.width = 0.1) +
    geom_parallel_sets_labels(color = "white")
cowplot::ggsave2(
  file.path(wd, "cluster_overlap_venn.pdf"),
  cluster_overlap_venn, width = 8, height = 10
)
