library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(egg)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("temporal_ordering", "tf_enrichment")
dir.create(wd, showWarnings = FALSE, recursive = TRUE)

# set directories, import files ------------------------------------------------
###############################################################################T

surface_fit <- syn("syn21444486") %>%
  read_csv()

temporal_ordering <- syn("syn21536903") %>%
  read_csv()

temporal_ordering_stats <- syn("syn21537190") %>%
  read_csv()

meta <- syn("syn21432975") %>%
  read_csv()

function_clusters <- syn("syn21478380") %>%
  read_csv()

condition_meta <- meta %>%
  distinct(condition, ERKi, Time, DOX)

pairwise_lfc <- syn("syn21432183") %>%
  read_csv()

lfc_long <- pairwise_lfc %>%
  # Remove genes that are NA in all conditions
  filter(
    !select(., starts_with("time")) %>%
      as.matrix() %>%
      is.na() %>%
      apply(1, all)
  ) %>%
  mutate_at(vars(starts_with("time")), replace_na, 0) %>%
  gather("condition", "log2FoldChange", -gene_name, -gene_id)


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
  summarize(targets = list(reduce(targets, union, .init = c()))) %>%
  ungroup()

# Load coexpression TF target data ---------------------------------------------
###############################################################################T

chea3_dbs <- tribble(
  ~tf_db, ~url,
  "gtex", "https://amp.pharm.mssm.edu/chea3/assets/tflibs/GTEx_Coexpression.gmt",
  "archs4", "https://amp.pharm.mssm.edu/chea3/assets/tflibs/ARCHS4_Coexpression.gmt",
  "lit_chip", "https://amp.pharm.mssm.edu/chea3/assets/tflibs/Literature_ChIP-seq.gmt"
) %>%
  mutate(
    tfs = map(
      url,
      ~.x %>%
        read_tsv(col_names = FALSE) %>%
        mutate(
          targets = pmap(., function(...) unname(na.omit(c(...)))[-1])
        ) %>%
        transmute(
          tf = str_replace(X1, fixed("_ARCHS4_PEARSON"), ""),
          targets
        )
    )
  ) %>%
  select(-url)

# Find combinations of TFs controlling bell targets ----------------------------
###############################################################################T

tf_dbs <- tribble(
  ~tf_db, ~tfs,
  "chea2", chea_tfs_agg
) %>%
  bind_rows(
    chea3_dbs
  )

# Using abs value of log2FoldChange for temporal induction ordering
temporal_ordering_abs <- temporal_ordering %>%
  filter(directed == "absolute")

asym_tfs <- tf_dbs %>%
  mutate(
    tfs = map(
      tfs,
      ~.x %>%
        inner_join(
          function_clusters %>%
            filter(class %in% c("high_erk", "low_erk")) %>%
            select(gene_name, tf_class = class, tf_direction = direction),
          by = c("tf" = "gene_name")
        ) %>%
        unnest_longer(targets)
    )
  )

all_tfs <- tf_dbs %>%
  mutate(
    tfs = map(
      tfs,
      ~.x %>%
        unnest_longer(targets)
    )
  )

bell_input_tfs <- function_clusters %>%
  filter(class == "bell") %>%
  inner_join(
    temporal_ordering_abs %>%
      filter(ERKi != "all") %>%
      select(gene_id, gene_name, ERKi, mid_induction, variance),
    by = c("gene_id", "gene_name")
  ) %>%
  inner_join(
    all_tfs %>%
      unnest(tfs),
    by = c("gene_name" = "targets")
  ) %>%
  inner_join(
    temporal_ordering_abs %>%
      filter(ERKi != "all") %>%
      select(gene_name, ERKi, mid_induction_tf = mid_induction, variance_tf = variance),
    by = c("tf" = "gene_name", "ERKi")
  )

bell_input_tfs_pass_agg <- bell_input_tfs %>%
  group_by(gene_id, gene_name, class, direction, class_combined, tf_db, tf) %>%
  summarize(
    mid_induction = weighted.mean(mid_induction, variance),
    mid_induction_tf = weighted.mean(mid_induction_tf, variance_tf),
    variance = max(variance),
    variance_tf = max(variance_tf)
  ) %>%
  ungroup() %>%
  filter(mid_induction_tf <= mid_induction) %>%
  # group_by(gene_name) %>%
  # filter(all(c("high_erk", "low_erk") %in% tf_class)) %>%
  # ungroup() %>%
  arrange(gene_name)

write_csv(
  bell_input_tfs_pass_agg,
  file.path(wd, "bell_genes_tfs_erki_aggregated.csv")
)

bell_input_tfs_pass <- bell_input_tfs %>%
  filter(mid_induction_tf <= mid_induction) %>%
  # group_by(gene_name, ERKi) %>%
  # filter(all(c("high_erk", "low_erk") %in% tf_class)) %>%
  # ungroup() %>%
  arrange(gene_name)

write_csv(
  bell_input_tfs_pass,
  file.path(wd, "bell_genes_tfs.csv")
)

