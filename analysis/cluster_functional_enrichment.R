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
library(Homo.sapiens)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("functional_enrichment")
dir.create(wd, recursive = TRUE, showWarnings = FALSE)

# set directories, import files ------------------------------------------------
###############################################################################T

function_clusters <- syn("syn21576614") %>%
  read_csv()

deseq_padj <- syn("syn21432184") %>%
  read_csv()

deseq_res <- syn("syn22686427") %>%
  read_csv()

cluster_names <- syn("syn21567677") %>%
  read_csv() %>%
  bind_rows(list(class_combined = "unknown", class_name = "Unknown")) %>%
  mutate_at(vars(class_name), . %>% as.factor() %>% fct_inorder())

temporal_ordering <- syn("syn21536903") %>%
  read_csv()

meta <- syn("syn21432975") %>%
  read_csv()

condition_meta <- meta %>%
  distinct(condition, DMSO, ERKi, Time, DOX)

surface_fit_p <- syn("syn22800020") %>%
  read_csv() %>%
  mutate(across(c(gene_id, gene_name), .fns = ~str_replace_all(.x, fixed("'"), "")))

# Metascape ----------------------------------------------------------------------
###############################################################################T

metascape_input <- function_clusters %>%
  select(gene_name, consensus) %>%
  filter(!consensus %in% c("no_response_0", "none")) %>%
  group_by(consensus) %>%
  mutate(id = 1:n()) %>%
  ungroup() %>%
  pivot_wider(id, names_from = consensus, values_from = gene_name, values_fill = "") %>%
  select(-id)

write_csv(metascape_input, file.path(wd, "metascape_input_consensus_sets.csv"))

# Metascape all conditions -----------------------------------------------------
###############################################################################T

metascape_input <- deseq_res %>%
  semi_join(
    condition_meta %>%
      filter(Time == 24),
    by = "condition"
  ) %>%
  filter(padj <= 0.05) %>%
  group_by(condition) %>%
  arrange(padj, .by_group = TRUE) %>%
  slice_head(n = 500) %>%
  distinct(gene_id, condition) %>%
  mutate(id = 1:n()) %>%
  ungroup() %>%
  pivot_wider(id, names_from = condition, values_from = gene_id, values_fill = "") %>%
  select(-id)

write_csv(metascape_input, file.path(wd, "metascape_input_all_conditions.csv"))

# Querying CMap ----------------------------------------------------------------
###############################################################################T

library(clueR)
library(biomaRt)

mart <- biomaRt::useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl"
)
ensembl_gene_id_mapping_biomart <- biomaRt::select(
  mart,
  unique(deseq_res$gene_id),
  c("entrezgene_id", "ensembl_gene_id"), "ensembl_gene_id"
) %>%
  as_tibble() %>%
  distinct()

clue_gmt <- deseq_res %>%
  # semi_join(
  #   condition_meta %>%
  #     filter(Time == 24 | (DOX == 1 & ERKi %in% c(0, 1000))),
  #   by = "condition"
  # ) %>%
  filter(padj <= 0.05) %>%
  arrange(padj, .by_group = TRUE) %>%
  inner_join(
    ensembl_gene_id_mapping_biomart,
    by = c("gene_id" = "ensembl_gene_id")
  ) %>%
  dplyr::transmute(
    gene_set = condition,
    gene_id = entrezgene_id,
    direction = if_else(log2FoldChange > 0, "up", "down")
  ) %>%
  mutate(
    chunk = gene_set %>%
      as.factor() %>%
      as.integer() %>%
      magrittr::divide_by_int(25)
  ) %>%
  split(.$chunk) %>%
  map(
    clueR::clue_gmt_from_df, drop_invalid = TRUE
  )

clue_jobs <- clue_gmt %>%
  imap(
    ~clueR::clue_query_submit(
      .x[["up"]], .x[["down"]],
      name = paste0("erk_", .y),
      use_fast_tool = FALSE
    )
  )

clue_result_files <- map_chr(
  clue_jobs, clue_query_download
)

clue_results <- clue_result_files %>%
  enframe("chunk", "result_path") %>%
  crossing(
    score_level = c("cell", "summary"),
    result_type = c("pert", "pcl")
  ) %>%
  mutate(
    data = pmap(
      .,
      function(result_path, score_level, result_type, ...)
        clue_parse_result(
          result_path, score_level = score_level, result_type = result_type,
          score_type = "tau"
        )
    )
  ) %>%
  dplyr::select(-chunk, -result_path) %>%
  unnest(data)

write_csv(
  clue_results,
  file.path(wd, "clue_results_all_conditions.csv.gz")
)

# Submit four clusters


clue_gmt <- tribble(
  ~direction, ~class_name, ~sign,
  "up", "bell", "-",
  "down", "bell", "+",
  "up", "full_range", "+",
  "down", "full_range", "-",
  "up", "high_erk", "+",
  "down", "high_erk", "-",
  "up", "low_erk", "-",
  "down", "low_erk", "+"
) %>%
  mutate(class = paste(class_name, sign, sep = "_")) %>%
  left_join(
    function_clusters %>%
      dplyr::select(gene_id, class = consensus),
    by = "class"
  ) %>%
  inner_join(
    surface_fit_p,
    by = "gene_id"
  ) %>%
  arrange(class_name, direction, fdr_lratio) %>%
  inner_join(
    ensembl_gene_id_mapping_biomart,
    by = c("gene_id" = "ensembl_gene_id")
  ) %>%
  dplyr::select(
    gene_set = class_name,
    gene_id = entrezgene_id,
    direction
  ) %>%
  mutate(
    chunk = gene_set %>%
      as.factor() %>%
      as.integer() %>%
      magrittr::divide_by_int(25)
  ) %>%
  split(.$chunk) %>%
  map(
    clueR::clue_gmt_from_df, drop_invalid = TRUE
  )

clue_jobs <- clue_gmt %>%
  imap(
    ~clueR::clue_query_submit(
      .x[["up"]], .x[["down"]],
      name = paste0("erk_", .y),
      use_fast_tool = FALSE
    )
  )

clue_result_files <- map_chr(
  clue_jobs, clue_query_download
)

clue_results <- clue_result_files %>%
  enframe("chunk", "result_path") %>%
  crossing(
    score_level = c("cell", "summary"),
    result_type = c("pert", "pcl")
  ) %>%
  mutate(
    data = pmap(
      .,
      function(result_path, score_level, result_type, ...)
        clue_parse_result(
          result_path, score_level = score_level, result_type = result_type,
          score_type = "tau"
        )
    )
  ) %>%
  dplyr::select(-chunk, -result_path) %>%
  unnest(data)

write_csv(
  clue_results,
  file.path(wd, "clue_results_consensus_clusters.csv.gz")
)


# Store to synapse -------------------------------------------------------------
###############################################################################T

activity <- Activity(
  "Connectivity map query",
  used = c(
    "syn21576614",
    "syn22686427",
    "syn22800020"
  ),
  executed = "https://github.com/clemenshug/erk_senescence/blob/master/analysis/cluster_functional_enrichment.R"
)

syn_cmap <- synExtra::synMkdir("syn21432134", "functional_enrichment", "cmap", .recursive = TRUE)

c(
  file.path(wd, "clue_results_consensus_clusters.csv.gz"),
  file.path(wd, "clue_results_all_conditions.csv.gz")
) %>%
  synStoreMany(syn_cmap, activity = activity)



# EnrichR ----------------------------------------------------------------------
###############################################################################T

function_clusters_enrichr <- function_clusters %>%
  select(-consensus_n) %>%
  gather("algorithm", "cluster", -gene_id, -gene_name) %>%
  group_by(cluster, algorithm) %>%
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

# Plotting functions -----------------------------------------------------------
###############################################################################T


plot_heatmap_gg <- function(df, aesthetic = aes(class, term, fill = signed_p), facet_by = NULL, ...) {
#  browser()
  facet_by_quo <- enquo(facet_by)
  abs_max <- df %>%
    pull(!!aesthetic[["fill"]]) %>%
    abs() %>%
    quantile(.95, names = FALSE, na.rm = TRUE) %>%
    round(1)
  mat <- df %>%
    dplyr::select(
      !!aesthetic[["x"]],
      !!aesthetic[["y"]],
      !!aesthetic[["fill"]],
      !!facet_by_quo
    ) %>%
    {mutate(., !!aesthetic[["x"]] := paste(!!aesthetic[["x"]], !!facet_by_quo)) %>% dplyr::select(-!!facet_by_quo)} %>%
    spread(!!aesthetic[["x"]], !!aesthetic[["fill"]], fill = 0) %>%
    column_to_rownames(quo_name(aesthetic[["y"]])) %>%
    as.matrix()
  row_clust <- hclust(dist(mat, method = "euclidian"), "ward.D2")
  df_ready <- df %>%
    mutate(
      !!aesthetic[["y"]] := factor(!!aesthetic[["y"]], levels = rownames(mat)[row_clust$order])
    )
  ggplot(df_ready, aesthetic) +
    geom_raster() +
    facet_wrap(vars(direction)) +
    scale_fill_distiller(palette = "RdBu", limits = c(-abs_max, abs_max), oob = scales::squish) +
    facet_wrap(vars(!!facet_by_quo)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1))
}


plot_heatmap_gg_bivariate <- function(
  df, fill_bi_var, aesthetic = aes(class, term, fill = p_signed), colormaps = NULL
) {
  fill_bi_var_quo <- enquo(fill_bi_var)
  fill_bi_vars <- df %>%
    pull(!!fill_bi_var_quo) %>%
    unique()
  abs_max <- df %>%
    pull(!!aesthetic[["fill"]]) %>%
    abs() %>%
    quantile(.95, names = FALSE, na.rm = TRUE) %>%
    round(1)
  mat <- df %>%
    dplyr::select(
      !!aesthetic[["x"]],
      !!aesthetic[["y"]],
      !!aesthetic[["fill"]],
      !!fill_bi_var_quo
    ) %>%
    mutate(!!aesthetic[["x"]] := paste(!!aesthetic[["x"]], !!fill_bi_var_quo)) %>%
    dplyr::select(-!!fill_bi_var_quo) %>%
    spread(!!aesthetic[["x"]], !!aesthetic[["fill"]]) %>%
    column_to_rownames(quo_name(aesthetic[["y"]])) %>%
    as.matrix()
  row_clust <- hclust(dist(mat, method = "euclidian"), "ward.D2")
  # browser()
  pals <- (if (!is.null(colormaps)) {
    colormaps
  } else {
    set_names(c("Reds", "Blues"), fill_bi_vars)
  }) %>%
    map(scales::col_numeric, domain = c(0, 1))
  # browser()
  df_ready <- df %>%
    dplyr::select(
      !!aesthetic[["x"]],
      !!aesthetic[["y"]],
      !!aesthetic[["fill"]],
      !!fill_bi_var_quo
    ) %>%
    mutate(
      !!aesthetic[["fill"]] := scales::rescale(!!aesthetic[["fill"]], to = c(0, 1), from = c(0, abs_max)) %>%
        magrittr::inset(. > 1, 1),
      !!aesthetic[["y"]] := factor(!!aesthetic[["y"]], levels = rownames(mat)[row_clust$order])
    ) %>%
    spread(!!fill_bi_var_quo, !!aesthetic[["fill"]]) %>%
    mutate(
      !!aesthetic[["fill"]] :=
        # colorspace::mixcolor(
        #   0.5,
        #   .[[fill_bi_vars[[1]]]] %>%
        #     pals[[fill_bi_vars[[1]]]]() %>%
        #     colorspace::hex2RGB(),
        #   .[[fill_bi_vars[[2]]]] %>%
        #     pals[[fill_bi_vars[[2]]]]() %>%
        #     colorspace::hex2RGB(),
        #   "LAB"
        # ) %>%
        new("RGB", coords = 1 - (1 - (
          .[[fill_bi_vars[[1]]]] %>%
            pals[[fill_bi_vars[[1]]]]() %>%
            colorspace::hex2RGB()
        )@coords +
          1 - (
            .[[fill_bi_vars[[2]]]] %>%
              pals[[fill_bi_vars[[2]]]]() %>%
              colorspace::hex2RGB()
          )@coords)) %>%
        colorspace::hex(fixup = TRUE)
    )
  ggplot(df_ready, aesthetic) +
    geom_raster() +
    scale_fill_identity()
}

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

# EnrichR heatmap --------------------------------------------------------------
###############################################################################T

clusters_go_enrichr_plot_data <- function_clusters_enrichr %>%
  unnest(enrichment) %>%
  group_nest(algorithm) %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        mutate(
          Term_combined = paste(str_sub(database, end = 4L), Term, sep = "_")
        ) %>%
        arrange(desc(Combined.Score)) %>%
        filter(
          Term_combined %in% (
            c(
              filter(., Adjusted.P.value <= 0.05) %>%
                group_by(cluster) %>%
                dplyr::slice(1:5) %>%
                ungroup() %>%
                pull(Term_combined),
              filter(., Adjusted.P.value <= 0.05) %>%
                pull(Term_combined)
            ) %>%
              unique() %>%
              head(50)
          )
        ) %>%
        mutate(
          neg_log10_p = -log10(Adjusted.P.value)
        ) %>%
        inner_join(
          cluster_names,
          by = c("cluster" = "class_combined")
        )
    )
  )

clusters_go_enrichr_plots <- clusters_go_enrichr_plot_data %>%
  mutate(
    data = map(
      data,
      ~plot_heatmap_gg(
        .x %>%
          mutate(Combined.Score = log2(Combined.Score)),
        aes(class_name, Term_combined, fill = Combined.Score), facet_by = NULL
      )
    )
  )

pwalk(
  clusters_go_enrichr_plots,
  function(algorithm, data, ...) {
    ggsave(
      file.path(wd, paste0("clusters_go_enrichr_heatmap_", algorithm, ".pdf")),
      data +
        ggtitle(algorithm),
      width = 10, height = 12
    )
  }
)

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

go_gene_mapping <- map(
  go_objects,
  ~usedGO(.x) %>%
    {genesInTerm(.x, .)} %>%
    enframe("id", "gene_id")
) %>%
  bind_rows() %>%
  unchop(gene_id) %>%
  genebabel::join_hgnc("gene_id", "ensembl_gene_id", c("symbol"))

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
  go_genes <- go_gene_mapping %>%
    filter(gene_id %in% gene_set) %>%
    group_by(id) %>%
    summarize(
      gene_symbols = paste(unique(symbol), collapse = "|"),
      gene_ids = paste(unique(gene_id), collapse = "|")
    ) %>%
    ungroup()
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
    dplyr::left_join(
      go_genes, by = "id"
    ) %>%
    dplyr::mutate_at(vars(pval), ~as.numeric(gsub("< 1e-30", "1e-30", .x)))
}

surface_fit_go <- topgo_enrichment(
  unique(function_clusters$gene_id), all_genes, "bp"
)

openxlsx::write.xlsx(
  surface_fit_go %>%
    arrange(pval),
  file.path(wd, "top_go_results_surface_fit_genes.xlsx")
)

# topGO on function clusters ---------------------------------------------------
###############################################################################T

plan(multisession(workers = 6))
function_clusters_topgo <- function_clusters %>%
  dplyr::select(-consensus_n) %>%
  gather("algorithm", "cluster", -gene_id, -gene_name) %>%
  filter(algorithm == "consensus") %>%
  group_nest(algorithm, cluster) %>%
  crossing(go_domain = c("bp", "mf")) %>%
  mutate(
    enrichment = future_map2(
      data, go_domain,
      ~topgo_enrichment(.x$gene_id, all_genes, go_domain = .y),
      .progress = TRUE
    )
  )

function_clusters_topgo_df <- function_clusters_topgo %>%
  dplyr::select(-data) %>%
  unnest(enrichment) %>%
  arrange(pval) %>%
  group_by(algorithm) %>%
  filter(
    id %in% (
      c(
        filter(., pval <= 0.05) %>%
          group_by(cluster) %>%
          dplyr::slice(1:5) %>%
          ungroup() %>%
          pull(id),
        filter(., pval <= 0.05) %>%
          pull(id)
      ) %>%
        unique() %>%
        head(50)
    )
  ) %>%
  ungroup()

write_csv(
  function_clusters_topgo_df,
  file.path(wd, paste("consensus_clusters_topgo_enrichment_table.csv"))
)

function_clusters_topgo_mat <- function_clusters_topgo_df %>%
  inner_join(cluster_names, by = c("cluster" = "class_combined")) %>%
  transmute(
    algorithm,
    cluster,
    term = paste0(go_domain, "_", term),
    pval = -log10(pval)
  ) %>%
  group_nest(algorithm) %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        spread(cluster, pval) %>%
        column_to_rownames("term") %>%
        as.matrix()
    )
  )

function_clusters_topgo_hm <- function_clusters_topgo_mat %>%
  mutate(
    data = map(data, plot_heatmap)
  )

pwalk(
  function_clusters_topgo_hm,
  function(algorithm, data) {
    ggsave(
      file.path(wd, paste0("go_heatmap_", algorithm, ".pdf")),
      data, width = 9, height = 12
    )
  }
)

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

plan(multisession(workers = 10))
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

x <- syn("syn21576614") %>%
  read_csv()

x %>% filter(
  pmap_lgl(
    list(function_fitting, k_medoids, linear_model),
    function(...) length(unique(list(...))) == 3
  )
)


temporal_ordering_abs_mid_ind_topgo_df <- temporal_ordering_abs_mid_ind_topgo %>%
  separate(mid_induction, c("mid_induction", "direction"), sep = "_") %>%
  filter(ERKi %in% c("0", "1000")) %>%
  dplyr::select(-data) %>%
  unnest(enrichment) %>%
  arrange(pval) %>%
  filter(
    id %in% (
      c(
        filter(., pval <= 0.05) %>%
          group_by(mid_induction, ERKi) %>%
          dplyr::slice(1:5) %>%
          ungroup() %>%
          pull(id),
        filter(., pval <= 0.05) %>%
          pull(id)
      ) %>%
        unique() %>%
        head(50)
    )
  ) %>%
  mutate(
    neg_log10_p = -log10(pval),
    signed_p = neg_log10_p * if_else(direction == "pos", 1, -1),
    term = paste0(go_domain, "_", term)
  ) %>%
  arrange(term, mid_induction, ERKi, direction)


temporal_ordering_abs_mid_ind_topgo_hm <- temporal_ordering_abs_mid_ind_topgo_df %>%
  mutate(
    mid_induction = factor(mid_induction, levels = c("early", "mid", "late"))
  ) %>%
  arrange(ERKi, mid_induction) %>%
  mutate(
    class = as.factor(
      paste("ERKi", ERKi, mid_induction, sep = "_")
    ) %>%
      fct_inorder()
  ) %>%
  plot_heatmap_gg_bivariate(
    fill_bi_var = direction,
    aesthetic = aes(class, term, fill = neg_log10_p),
    colormaps = list(
      # "pos" = "Reds", "neg" = "Blues"
      "pos" = c("#FFFFFFFF", "#DD000000"),
      "neg" = c("#FFFFFFFF", "#0000DD00")
    )
  ) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(fill = "signed\n-log10(p)")
ggsave(
  file.path(wd, "temporal_ordering_mid_induction_go_heatmap_bivar.pdf"),
  temporal_ordering_abs_mid_ind_topgo_hm, width = 10, height = 8
)



# topGO on high/low ERK responders by time -------------------------------------
###############################################################################T


cons_clusters <- function_clusters %>%
  dplyr::select(
    gene_id, gene_name, consensus_n, consensus
  ) %>%
  extract(
    consensus, c("cluster", "direction"),
    regex = "^(full_range|high_erk|low_erk|no_response|bell)_([0\\+-])$", remove = FALSE
  )

high_low_clusters <- cons_clusters %>%
  filter(cluster %in% c("high_erk", "low_erk"))

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

plan(multisession(workers = 10))
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


high_low_temp_clusters <- temporal_ordering %>%
  filter(directed == "directed", ERKi %in% c(0, 1000)) %>%
  dplyr::select(gene_id, mid_induction, ERKi) %>%
  mutate(high_low = if_else(ERKi == 0, "low_erk", "high_erk")) %>%
  inner_join(high_low_clusters, by = c("gene_id", "high_low" = "cluster"))

plan(multisession(workers = 4))
high_low_clusters_go <- high_low_temp_clusters %>%
  group_nest(high_low, mid_induction, direction) %>%
  mutate_at(vars(mid_induction), as.character) %>%
  bind_rows(
    mutate(
      .,
        mid_induction = cut(as.integer(mid_induction), breaks = c(0, 3, 10, Inf), labels = c("early", "mid", "late"))
      ) %>%
      group_by(high_low, mid_induction, direction) %>%
      summarize(data = list(bind_rows(data))) %>%
      ungroup()
  ) %>%
  crossing(go_domain = c("bp", "mf")) %>%
  mutate(
    enrichment = future_map2(
      data, go_domain,
      ~topgo_enrichment(.x$gene_id, all_genes, go_domain = .y),
      .progress = TRUE
    )
  )

write_rds(
  high_low_clusters_go,
  file.path(wd, "high_low_erk_temporal_go_enrichment.rds"),
  compress = "gz"
)

# high_low_clusters_go <- read_rds(file.path(wd, "high_low_erk_temporal_go_enrichment.rds"))


high_low_clusters_go_plot_data <- high_low_clusters_go %>%
  mutate(
    direction = if_else(xor(direction == "+", high_low == "low_erk"), "up", "down")
  ) %>%
  filter(mid_induction %in% c("early", "mid", "late")) %>%
  dplyr::select(-data) %>%
  unnest(enrichment) %>%
  arrange(pval) %>%
  group_nest(high_low) %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        filter(
          id %in% (
            c(
              filter(., pval <= 0.05) %>%
                group_by(mid_induction, direction) %>%
                dplyr::slice(1:5) %>%
                ungroup() %>%
                pull(id),
              filter(., pval <= 0.05) %>%
                pull(id)
            ) %>%
              unique() %>%
              head(50)
          )
        ) %>%
        mutate(
          neg_log10_p = -log10(pval),
          signed_p = neg_log10_p * if_else(direction == "up", 1, -1),
          term = paste0(go_domain, "_", term),
          mid_induction = factor(mid_induction, levels = c("early", "mid", "late"))
        ) %>%
        arrange(term, mid_induction, direction)
    )
  )

high_low_clusters_go_plot_data %>%
  unnest(data) %>%
  arrange(pval) %>%
  write_csv(file.path(wd, "high_low_erk_temporal_go_enrichment_top.csv"))

high_low_clusters_go_plot <- high_low_clusters_go_plot_data %>%
  mutate(
    data = map(
      data,
      ~plot_heatmap_gg(
        .x,
        aes(mid_induction, term, fill = signed_p), facet_by = direction
      )
    )
  )

pwalk(
  high_low_clusters_go_plot,
  function(high_low, data, ...) {
    ggsave(
      file.path(wd, paste0("high_low_temporal_ordering_go_heatmap", high_low, ".pdf")),
      data +
        ggtitle(high_low)
    )
  }
)

high_low_clusters_go_enrichr <- high_low_temp_clusters %>%
  group_nest(high_low, mid_induction, direction) %>%
  mutate_at(vars(mid_induction), as.character) %>%
  bind_rows(
    mutate(
      .,
      mid_induction = cut(as.integer(mid_induction), breaks = c(0, 3, 10, Inf), labels = c("early", "mid", "late"))
    ) %>%
      group_by(high_low, mid_induction, direction) %>%
      summarize(data = list(bind_rows(data))) %>%
      ungroup()
  ) %>%
  mutate(
    enrichment = map(
      data,
      ~enrichr(
        .x$gene_name,
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
  )

high_low_clusters_go_enrichr_plot_data <- high_low_clusters_go_enrichr %>%
  mutate(
    direction = if_else(xor(direction == "+", high_low == "low_erk"), "up", "down")
  ) %>%
  filter(mid_induction %in% c("early", "mid", "late")) %>%
  dplyr::select(-data) %>%
  unnest(enrichment) %>%
  mutate(
    Term_combined = paste(str_sub(database, end = 4L), Term, sep = "_")
  ) %>%
  arrange(desc(Combined.Score)) %>%
  group_nest(high_low) %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        filter(
          Term_combined %in% (
            c(
              filter(., Adjusted.P.value <= 0.05) %>%
                group_by(mid_induction, direction) %>%
                dplyr::slice(1:5) %>%
                ungroup() %>%
                pull(Term_combined),
              filter(., Adjusted.P.value <= 0.05) %>%
                pull(Term_combined)
            ) %>%
              unique() %>%
              head(50)
          )
        ) %>%
        mutate(
          neg_log10_p = -log10(Adjusted.P.value),
          signed_p = neg_log10_p * if_else(direction == "up", 1, -1),
          mid_induction = factor(mid_induction, levels = c("early", "mid", "late"))
        )
    )
  )

high_low_clusters_go_enrichr_plot <- high_low_clusters_go_enrichr_plot_data %>%
  mutate(
    data = map(
      data,
      ~plot_heatmap_gg(
        .x,
        aes(mid_induction, Term_combined, fill = signed_p), facet_by = direction
      )
    )
  )

pwalk(
  high_low_clusters_go_enrichr_plot,
  function(high_low, data, ...) {
    ggsave(
      file.path(wd, paste0("high_low_temporal_ordering_go_enrichr_heatmap", high_low, ".pdf")),
      data +
        ggtitle(high_low)
    )
  }
)

