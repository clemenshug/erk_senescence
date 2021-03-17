library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(egg)
library(seriation)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("temporal_ordering")
dir.create(wd, showWarnings = FALSE)

# set directories, import files ------------------------------------------------
###############################################################################T

function_clusters <- syn("syn21576614") %>%
  read_csv()
deseq_res <- syn("syn21432187") %>%
  read_rds()
pairwise_lfc <- syn("syn21432183") %>%
  read_csv()
pairwise_padj <- syn("syn21432184") %>%
  read_csv()
meta <- syn("syn21432975") %>%
  read_csv()

condition_meta <- meta %>%
  distinct(condition, DMSO, ERKi,  Time, DOX)


# CLuster high/low ERK classes -------------------------------------------------
###############################################################################T

lfc_long <- pairwise_lfc %>%
  gather("condition", "lfc", -starts_with("gene"))

cons_clusters <- function_clusters %>%
  select(
    gene_id, gene_name, consensus_n, consensus
  ) %>%
  extract(
    consensus, c("cluster", "direction"),
    regex = "^(full_range|high_erk|low_erk|no_response|bell)_([0\\+-])$", remove = FALSE
  )

high_low_clusters <- cons_clusters %>%
  filter(cluster %in% c("high_erk", "low_erk"))

high_low_conditions <- condition_meta %>%
  mutate(
    high_low = case_when(
      ERKi == 0 & DOX == 1 ~ "high_erk",
      ERKi == 1000 & DOX == 1 ~ "low_erk",
      TRUE ~ "neither"
    )
  ) %>%
  filter(high_low != "neither") %>%
  group_nest(high_low) %>%
  with(set_names(data, high_low))

set.seed(42)
high_low_lfc <- high_low_clusters %>%
  group_nest(consensus, cluster, direction) %>%
  mutate(
    data = map2(
      data, cluster,
      ~high_low_conditions[[.y]] %>%
        select(condition, Time) %>%
        inner_join(
          lfc_long, by = "condition"
        ) %>%
        filter(gene_id %in% .x[["gene_id"]]) %>%
        arrange(Time) %>%
        mutate(Time = fct_inorder(as.character(Time))) %>%
        select(Time, gene_id, lfc) %>%
        spread(Time, lfc, fill = 0) %>%
        mutate_at(vars(-gene_id), replace_na, 0) %>%
        column_to_rownames("gene_id") %>%
        as.matrix()
    )
  )

high_low_hc <- high_low_lfc %>%
  mutate(
    clust = map(
      data,
      function(m) {
        d <- cor(t(m), method = "kendall") %>%
          {1 - .} %>%
          as.dist()
        # browser()
        hclust(d, method = "average")
          # reorder(dist = d)
      }
    )
  )

high_low_hc_plots <- high_low_hc %>%
  mutate(
    data = map2(
      clust, data,
      function(h, m) {
        x <- m %>%
          as_tibble(rownames = "gene_id")
        # browser()
        x %>%
          mutate(gene_id = factor(gene_id, levels = h$labels)) %>%
          gather("time", "lfc", -gene_id) %>%
          # group_by(gene_id) %>%
          # mutate(lfc = scale(lfc)) %>%
          # ungroup() %>%
          mutate(time = factor(time, levels = as.character(sort(unique(as.integer(time)))))) %>%
          ggplot(aes(time, gene_id, fill = lfc)) +
            geom_raster() +
            scale_fill_distiller(palette = "RdBu", limits = c(-3, 3), oob = scales::squish)
      }
    )
  )

set.seed(42)
high_low_kmed <- high_low_lfc %>%
  mutate(
    clust = map(
      data,
      function(m) {
        d <- cor(t(m), method = "kendall") %>%
          {1 - .} %>%
          as.dist()
        cluster::pam(d, k = 4)
      }
    ),
    clust_df = map(
      clust,
      ~.x %>%
        chuck("clustering") %>%
        enframe("gene_id", "cluster")
    )
  )


high_low_kmed_plots <- high_low_kmed %>%
  mutate(
    data = map2(
      clust_df, data,
      function(h, m) {
        m %>%
          as_tibble(rownames = "gene_id") %>%
          gather("time", "lfc", -gene_id) %>%
          group_by(gene_id) %>%
          mutate(lfc = scale(lfc)) %>%
          ungroup() %>%
          inner_join(h, by = "gene_id") %>%
          mutate(time = factor(time, levels = as.character(sort(unique(as.integer(time)))))) %>%
          ggplot(aes(time, lfc, group = gene_id)) +
            geom_line() +
            facet_wrap(~cluster)
      }
    )
  )
