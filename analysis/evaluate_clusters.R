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

syn_clustering <- "syn21432135"

wd <- here("clustering", "overlap")
dir.create(wd, recursive = TRUE, showWarnings = FALSE)

# set directories, import files ------------------------------------------------
###############################################################################T

km_clusters <- syn("syn21574266") %>%
  read_csv()
lm_clusters <- syn("syn21448572") %>%
  read_csv() %>%
  rename(cluster = class, cluster_combined = class_combined) %>%
  select(-low, -high)
function_clusters <- syn("syn21478380") %>%
  read_csv() %>%
  rename(cluster = class, cluster_combined = class_combined)
surface_fit <- syn("syn21444486") %>%
  read_csv()

# Merge clustering data --------------------------------------------------------
###############################################################################T

cluster_algos <- tribble(
  ~algorithm, ~clusters,
  "k_medoids", km_clusters,
  "linear_model", lm_clusters %>%
    filter(gene_id %in% surface_fit$gene_id),
  "function_fitting", function_clusters
)

cluster_overlap <- cluster_algos %>%
  unnest(clusters) %>%
  select(algorithm, gene_id, cluster_combined) %>%
  inner_join(
    surface_fit %>%
      distinct(gene_id, gene_name),
    by = "gene_id"
  ) %>% {
    inner_join(
      spread(., algorithm, cluster_combined) %>%
        drop_na(),
      group_by(., gene_name, gene_id) %>%
        summarize(
          ctable = list(table(cluster_combined)),
        ) %>%
        ungroup() %>%
        transmute(
          gene_name, gene_id,
          consensus = map(
            ctable,
            function(t) {
              h <- order(t, decreasing = TRUE)[[1]]
              if (t[[h]] < 2)
                list(consensus = "none", consensus_n = NA_integer_)
              else
                list(consensus = names(t)[[h]], consensus_n = t[[h]])
            }
          )
        )
    )
  } %>%
  unnest_wider(consensus)

write_csv(
  cluster_overlap,
  file.path(wd, "cluster_overlap.csv")
)

cluster_overlap_venn <- cluster_overlap %>%
  group_by_at(vars(3:6)) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  gather_set_data(1:4) %>%
  ggplot(aes(x, id = id, split = y, value = n)) +
    geom_parallel_sets(aes(fill = k_medoids), alpha = 0.3, axis.width = 0.1) +
    geom_parallel_sets_axes(axis.width = 0.1) +
    geom_parallel_sets_labels(color = "white")
cowplot::ggsave2(
  file.path(wd, "cluster_overlap_venn.pdf"),
  cluster_overlap_venn, width = 8, height = 12
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

activity <- Activity(
  "Evaluate and compare different clustering methods",
  used = c(
    "syn21574266",
    "syn21448572",
    "syn21478380",
    "syn21444486"
  ),
  executed = "https://github.com/clemenshug/erk_senescence/blob/master/analysis/evaluate_clusters.R"
)

syn_overlap <- Folder("overlap", parent = syn_clustering) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(wd, "cluster_overlap.csv")
) %>%
  synStoreMany(syn_overlap, activity = activity)

