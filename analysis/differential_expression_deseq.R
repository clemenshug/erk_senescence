library(tidyverse)
library(DESeq2)
library(readxl)
library(furrr)
library(rtracklayer)

synapser::synLogin()
syn <- synExtra::synDownloader("~/data", .cache = TRUE)


meta <- syn("syn21411559") %>%
  read_excel() %>%
  mutate(
    condition = paste0("time",Time,"dox",DOX) %>%
      as.factor() %>%
      relevel("time0dox0"),
    ERKiConc = as.factor(ERKi) %>%
      relevel("0"),
    Repeat = as.factor(Repeat),
    condition_all = paste0("time",Time,"dox",DOX,"conc",ERKi) %>%
      as.factor() %>%
      relevel("time0dox0conc0"),
    Sample_ID = paste0("S", as.integer(str_extract_all(Sample_ID, "[0-9]+")))
  )

gene_map <- syn("syn21411554") %>%
  read_csv()
gene_map_protein_coding <- gene_map %>%
  filter(gene_biotype == "protein_coding")

countData <- syn("syn21411551") %>%
  read_tsv() %>%
  filter(id %in% gene_map_protein_coding$gene_id)

dds <- DESeqDataSetFromMatrix(
  countData = countData %>%
    column_to_rownames("id") %>%
    as.matrix() %>%
    .[, meta$Sample_ID],
  colData = meta %>%
      column_to_rownames("Sample_ID"),
  design = ~ condition_all + Repeat
)
dds

de <- DESeq(dds)

dir.create("deseq_protein_coding", showWarnings = FALSE)

write_rds(de, file.path("deseq_protein_coding", "deseq_pairwise.rds"))

resultsNames(de)

# Compare specifically two conditions of interest for QC of the ERKi
# qc_res <- lfcShrink(de, contrast = c("condition_all", "time24dox1conc1000", "time24dox0conc250"), type = "ashr")
qc_res <- results(de, contrast = c("condition_all", "time24dox1conc1000", "time24dox0conc250"))
write_csv(
  qc_res %>%
    as.data.frame() %>%
    rownames_to_column("gene_id") %>%
    left_join(select(gene_map, gene_id, gene_name)),
  file.path("deseq", paste0("deseq_erki_qc_comparison_time24dox1conc1000_vs_time24dox0conc250", ".csv.gz"))
)


plan(multisession(workers = 4))
all_results <- resultsNames(de) %>%
  keep(str_starts, "condition") %>%
  set_names() %>%
  future_map(
    ~lfcShrink(de, coef = .x, type = "apeglm"),
    .progress = TRUE
  ) %>%
  map(as.data.frame) %>%
  map(rownames_to_column, "gene_id")

all_results_table <- all_results %>%
  enframe("condition", "res") %>%
  mutate(
    condition = str_match(condition, "condition_all_(.+)_vs_time0dox0conc0")[, 2]
  ) %>%
  unnest(res)

write_csv(
  all_results_table %>%
    left_join(select(gene_map, gene_id, gene_name)),
  file.path("deseq", paste0("deseq_pairwise_results", ".csv.gz"))
)

synStoreMany(
  file.path("deseq", paste0("deseq_pairwise_results", ".csv.gz")), "syn21432177"
)

all_results_tables
  select(condition, gene_id, log2FoldChange, lfcSE, padj) %>%
  gather("variable", "value", log2FoldChange, lfcSE, padj) %>%
  group_nest(variable) %>%
  mutate(
    data = map(
      data,
      ~spread(.x, condition, value) %>%
        left_join(select(gene_map, gene_id, gene_name)) %>%
        select(gene_id, gene_name, everything())
    )
  )

pwalk(
  all_results_tables,
  function(variable, data, ...)
    write_csv(
      data,
      file.path("deseq_protein_coding", paste0("deseq_pairwise_", variable, ".csv"))
    )
)


# deseq_datasets <- meta %>%
#   group_nest(ERKi, .key="meta_table") %>%
#   mutate(
#     meta_table = map(meta_table, bind_rows, meta %>% filter(condition == "time0dox0")) %>%
#       map(distinct)
#   ) %>%
#   filter(ERKi != 0) %>%
#   mutate(
#     count_table = map(meta_table, ~countData[, rownames(.x)])
#   )

plan(multisession(workers = 6))
deseq_linear <- meta %>%
  mutate(log2ERKi = log2(ERKi + 1)) %>%
  group_nest(condition, .key= "meta") %>%
  filter(condition != "time0dox0") %>%
  mutate(count = map(meta, ~select(countData, id, one_of(.x$Sample_ID)))) %>%
  arrange(condition) %>%
  mutate(
    deseq = future_map2(
      meta, count,
      ~DESeqDataSetFromMatrix(
        countData = .y %>%
          column_to_rownames("id") %>%
          as.matrix() %>%
          .[, .x$Sample_ID],
        colData = column_to_rownames(.x, "Sample_ID"),
        design = ~ log2ERKi + Repeat) %>% DESeq(), .progress = TRUE
      )
  )

deseq_linear_res <- deseq_linear %>%
  mutate(res = map(deseq, results, name = "log2ERKi") %>% map(as.data.frame) %>% map(rownames_to_column, "gene_id"))

deseq_linear_res_tables <- deseq_linear_res %>%
  select(condition, res) %>%
  unnest(res) %>%
  select(condition, gene_id, log2FoldChange, lfcSE, padj) %>%
  gather("variable", "value", log2FoldChange, lfcSE, padj) %>%
  group_nest(variable) %>%
  mutate(
    data = map(
      data,
      ~spread(.x, condition, value)
    )
  )

pwalk(
  deseq_linear_res_tables,
  function(variable, data, ...)
    write_csv(
      data %>%
        left_join(select(gene_map, gene_id, gene_name)) %>%
        select(gene_id, gene_name, everything()),
      file.path("deseq_protein_coding", paste0("deseq_linear_", variable, ".csv"))
    )
)

pwalk(
  deseq_linear_res,
  function(condition, res, ...)
    write_csv(
      res %>%
        left_join(select(gene_map, gene_id, gene_name)) %>%
        select(gene_id, gene_name, everything()),
      file.path("deseq_protein_coding", paste0("deseq_linear_", condition, ".csv"))
    )
)
