library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(RColorBrewer)
library(pheatmap)
library(ggforce)
library(viridis)
library(fgsea)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("gsea_enrichment")
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

# Gene sets --------------------------------------------------------------------
###############################################################################T

library(msigdbr)

of_interest <- list(
  "REACTOME_ONCOGENE_INDUCED_SENESCENCE",
  "BIOCARTA_ERK_PATHWAY",
  "REACTOME_ERK_MAPK_TARGETS",
  "REACTOME_SIGNALLING_TO_ERKS",
  "ST_ERK1_ERK2_MAPK_PATHWAY",
  "REACTOME_CELLULAR_SENESCENCE",
  "REACTOME_MAPK_FAMILY_SIGNALING_CASCADES",
  "REACTOME_MAPK1_ERK2_ACTIVATION",
  "REACTOME_MAPK3_ERK1_ACTIVATION",
  "HP_MELANOMA"
) %>%
  rlang::set_names() %>%
  c(
    list(
      "MAPK activation combined" = c(
        "REACTOME_MAPK1_ERK2_ACTIVATION",
        "REACTOME_MAPK3_ERK1_ACTIVATION"
      ),
      "MAPK combined" = c(
        "BIOCARTA_ERK_PATHWAY",
        "REACTOME_ERK_MAPK_TARGETS",
        "REACTOME_SIGNALLING_TO_ERKS",
        "ST_ERK1_ERK2_MAPK_PATHWAY",
        "REACTOME_MAPK1_ERK2_ACTIVATION",
        "REACTOME_MAPK3_ERK1_ACTIVATION"
      )
    )
  )

all_msigdbr <- msigdbr()

braf_oe_geo <- map(
  c("up", "down") %>%
    set_names(paste0("braf_v600e_geo_", .)),
  function(direction) {
    read_tsv(
      here(
        "functional_enrichment",
        paste0(
          "BRAF_overexpression%20human%20GSE46801%20sample%201951_from_Gene_Perturbations_from_GEO_", direction, ".gmt"
        )
      ),
      col_names = FALSE
    ) %>%
      {
        tibble(
          gene_symbol = as.character(.[1, ])[-c(1, 2)]
        )
      }
  }
)

msigdbr_of_interest <- of_interest %>%
  map(
    ~all_msigdbr %>%
      filter(gs_name %in% .x) %>%
      distinct(entrez_gene, gene_symbol = human_gene_symbol)
  ) %>%
  c(braf_oe_geo) %>%
  enframe("gene_set_name", "data")



parallel_param <- BiocParallel::MulticoreParam(4)

msigdbr_fgsea_res <- fgseaSimple(
  msigdbr_of_interest %>%
    {
      set_names(
        .[["data"]] %>%
          map("gene_symbol"),
        .[["gene_set_name"]]
      )
    },
  surface_fit_p %>%
    group_by(gene_name) %>%
    summarize(fdr_lratio = exp(mean(log(fdr_lratio))), .groups = "drop") %>%
    drop_na() %>%
    {
      set_names(
        .[["fdr_lratio"]] %>%
          pmin(1) %>%
          pmax(min(.[. > 0]) / 2),
          # log2() %>%
          # magrittr::multiply_by(-1),
        .[["gene_name"]]
      )
    },
  nperm = 1000,
  BPPARAM = parallel_param,
  scoreType = "pos"
)

fisher_enrich <- function(
  significant_symbols, gene_set_symbols, gene_universe
){
  fisher_table <- table(
    significant_gene = gene_universe %in% significant_symbols,
    gene_set_gene = gene_universe %in% gene_set_symbols
  )
  fisher_res <- fisher_table %>%
    fisher.test(
      alternative = "greater",
      conf.int = TRUE,
      conf.level = 0.95
    ) %>%
    broom::tidy()
  list(
    fisher_table = fisher_table,
    fisher_res = fisher_res
  )
}

library(magrittr)
fisher_enrich_bg <- function(
  n_significant_symbols, gene_set_symbols, gene_universe, reps = 1000
){
  significant_vec <- rep_len(FALSE, length(gene_universe))
  significant_vec[seq_len(n_significant_symbols)] <- TRUE
  replicate(
    reps,
    fisher_table <- table(
      significant_gene = sample(significant_vec),
      gene_set_gene = gene_universe %in% gene_set_symbols
    ) %>%
      fisher.test(
        alternative = "greater",
        conf.int = TRUE,
        conf.level = 0.95
      ) %>%
      broom::tidy(),
    simplify = FALSE
  ) %>%
    bind_rows()
}

surface_p_clean <- surface_fit_p %>%
  group_by(gene_name) %>%
  summarize(fdr_lratio = exp(mean(log(fdr_lratio))), .groups = "drop") %>%
  drop_na() %>%
  mutate(
    significant = fdr_lratio < 1e-20
  )

surface_p_clean_sig <- surface_p_clean %>%
  filter(significant)

genesets_of_interest <- function_clusters %>%
  group_nest(consensus, .key = "query_gene_set") %>%
  mutate(across(query_gene_set, as.list)) %>%
  inner_join(cluster_names, by = c("consensus" = "class_combined")) %>%
  rename(query_gene_set_name = class_name) %>%
  bind_rows(
    tibble(
      query_gene_set_name = "all significant",
      query_gene_set = list(surface_p_clean_sig)
    )
  )

library(furrr)

plan(multisession(workers = 8))

msigdbr_fisher_res <- msigdbr_of_interest %>%
  rename(gene_set = data) %>%
  crossing(genesets_of_interest) %>%
  {
    res <- map2(
      .[["query_gene_set"]], .[["gene_set"]],
      ~fisher_enrich(
        .x[["gene_name"]],
        .y[["gene_symbol"]],
        surface_p_clean[["gene_name"]]
      )
    )
    bg <- future_map2(
      .[["query_gene_set"]], .[["gene_set"]],
      ~fisher_enrich_bg(
        nrow(.x),
        .y[["gene_symbol"]],
        surface_p_clean[["gene_name"]],
        reps = 1000
      ),
      .progress = TRUE,
      .options = furrr_options(seed = TRUE)
    )
    # bg <- map2(
    #   .[["query_gene_set"]], .[["gene_set"]],
    #   ~fisher_enrich_bg(
    #     nrow(.x),
    #     .y[["gene_symbol"]],
    #     surface_p_clean[["gene_name"]],
    #     reps = 5
    #   )
    # )
    mutate(
      .,
      fisher_table = map(res, "fisher_table"),
      fisher_res = map(res, "fisher_res"),
      fisher_bg = map(bg, "p.value")
    )
  } %>%
  bind_cols(
    .[["fisher_res"]] %>%
      bind_rows()
  ) %>%
  arrange(p.value)


msigdbr_fisher_plot_dir <- file.path(wd, "msigdbr_fisher_plots")
dir.create(msigdbr_fisher_plot_dir)

msigdbr_fisher_res_table <- msigdbr_fisher_res %>%
  filter(!query_gene_set_name %in% c("other")) %>%
  head(n = 20) %>%
  select(
    query_gene_set_name, gene_set_name,
    p.value, odds_ratio = estimate
  ) %>%
  tableGrob()

ggsave(
  file.path(msigdbr_fisher_plot_dir, "msigdbr_fisher_res_table.pdf"),
  msigdbr_fisher_res_table,
  width = 10, height = 7
)

msigdbr_fisher_res_tables <- msigdbr_fisher_res %>%
  group_nest(query_gene_set_name, keep = TRUE) %>%
  filter(!query_gene_set_name %in% c("other")) %>%
  rowwise() %>%
  mutate(
    table = data %>%
      head(n = 20) %>%
      select(
        query_gene_set_name, gene_set_name,
        p.value, odds_ratio = estimate
      ) %>%
      tableGrob() %>%
      list()
  )

pwalk(
  msigdbr_fisher_res_tables,
  function(table, query_gene_set_name, ...) {
    ggsave(
      file.path(msigdbr_fisher_plot_dir, paste0("msigdbr_fisher_res_table_", query_gene_set_name, ".pdf")),
                table,
                width = 10, height = 7
      )
  }
)

msigdbr_fisher_res_tables <- msigdbr_fisher_res %>%
  group_nest(gene_set_name, keep = TRUE) %>%
  rowwise() %>%
  mutate(
    table = data %>%
      filter(!query_gene_set_name %in% c("other")) %>%
      head(n = 20) %>%
      select(
        query_gene_set_name, gene_set_name,
        p.value, odds_ratio = estimate
      ) %>%
      tableGrob() %>%
      list()
  )

pwalk(
  msigdbr_fisher_res_tables,
  function(table, gene_set_name, ...) {
    ggsave(
      file.path(msigdbr_fisher_plot_dir, paste0("msigdbr_fisher_res_table_", gene_set_name, ".pdf")),
      table,
      width = 10, height = 7
    )
  }
)

x <- msigdbr_fisher_res %>%
  arrange(p.value) %>%
  # head(n = 5) %>%
  pmap(
    function(p.value, gene_set_name, query_gene_set_name, fisher_table, fisher_bg, ...) {
      dist_plot <- ggplot(
        tibble(p.value = fisher_bg),
        aes(x = p.value)
      ) +
        geom_density(aes(y = stat(scaled))) +
        geom_segment(
          aes(x = x, xend = x, y = 0, yend = 1),
          data = tibble(x = p.value),
          color = "maroon"
        ) +
        scale_x_reverse() +
        labs(
          x = "p value",
          y = "density estimate"
        )

    }
  )



  fgseaSimple(
  msigdbr_of_interest %>%
    {
      set_names(
        .[["data"]] %>%
          map("gene_symbol"),
        .[["gene_set_name"]]
      )
    },
  surface_fit_p %>%
    group_by(gene_name) %>%
    summarize(fdr_lratio = exp(mean(log(fdr_lratio))), .groups = "drop") %>%
    drop_na() %>%
    {
      set_names(
        .[["fdr_lratio"]] %>%
          pmin(1) %>%
          pmax(min(.[. > 0]) / 2),
        # log2() %>%
        # magrittr::multiply_by(-1),
        .[["gene_name"]]
      )
    },
  nperm = 1000,
  BPPARAM = parallel_param,
  scoreType = "pos"
)
