library(tidyverse)
library(DESeq2)
library(readxl)
library(furrr)
library(rtracklayer)
library(here)
library(synapser)
library(synExtra)

wd <- here("deseq")
dir.create(wd, showWarnings = TRUE)

syn <- synDownloader(here("tempdl"), followLink = TRUE)

# Loading data -----------------------------------------------------------------
###############################################################################T

de <- syn("syn21432187") %>%
  read_rds()

# Normalize counts -------------------------------------------------------------
###############################################################################T

norm_counts <- counts(de, normalized = TRUE)
varstab <- varianceStabilizingTransformation(de, blind = FALSE) %>%
  assay()

combined <- tribble(
  ~type, ~data,
  "normalized", norm_counts,
  "variance_stabilized", varstab
) %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        as_tibble(rownames = "ensembl_gene_id") %>%
        gather("sample_id", "count", -ensembl_gene_id)
    )
  ) %>%
  unnest(data)

write_csv(
  combined,
  file.path(wd, "normalized_counts.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

activity <- Activity(
  "Normalize counts",
  used = c(
    "syn21432187"
  ),
  executed = "https://github.com/clemenshug/erk_senescence/blob/master/analysis/normalize_counts.R"
)

c(
  file.path(wd, "normalized_counts.csv.gz")
) %>%
  synStoreMany("syn21432177", activity = activity)

