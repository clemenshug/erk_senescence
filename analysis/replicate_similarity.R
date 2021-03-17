library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(DESeq2)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("replicate_analysis")
dir.create(wd, recursive = TRUE, showWarnings = FALSE)

# set directories, import files ------------------------------------------------
###############################################################################T

raw_counts <- syn("syn21411551") %>%
  read_tsv()

meta <- syn("syn21432975") %>%
  read_csv() %>%
  group_by(condition) %>%
  mutate(
    replicate = seq_len(n())
  ) %>%
  ungroup()

deseq_pairwise <- syn("syn21432187") %>%
  read_rds()

unloadNamespace("synExtra")
unloadNamespace("synapser")
unloadNamespace("PythonEmbedInR")

# Normalize counts -------------------------------------------------------------
###############################################################################T

counts_vst <- deseq_pairwise %>%
  varianceStabilizingTransformation() %>%
  assay()

# PCA with normalized counts ---------------------------------------------------
###############################################################################T

pca_data <- counts_vst %>%
  prcomp()


pca_plot <- function (data, meta, aes = ggplot2::aes(PC1, PC2), extra_layers = NULL, ...) {
  p <- prcomp(data, ...)
  pstats <- t(summary(p)$importance) %>%
    tibble::as_tibble(rownames = "component") %>%
    dplyr::rename(sdev = `Standard deviation`, prop_var = `Proportion of Variance`, cum_prop = `Cumulative Proportion`) %>%
    dplyr::mutate(component = factor(component, levels = unique(component))) %>%
    # Remove all components after cumsum reaches .999
    dplyr::filter(cumsum(dplyr::lag(cum_prop > .95, default = FALSE)) <= 1) %>%
    # Maximum 10 components
    dplyr::slice(1:min(10, n()))
  ploadings <- p$x %>%
    as.data.frame() %>%
    tibble::rownames_to_column("condition") %>%
    tibble::as_tibble() %>%
    dplyr::inner_join(meta, by = "condition")
  p_plot <- ggplot(ploadings, aes)
  if (!is.null(extra_layers))
    p_plot <- p_plot + extra_layers
  # var_plot <- ggplot(pstats, ggplot2::aes(x = 1, y = prop_var, fill = component)) +
  #   geom_col(position = "stack") +
  #   # geom_text(ggplot2::aes(y = cum_prop, label = prop_var), halign = 0, valign = 1) +
  #   coord_flip() +
  #   guides(fill = FALSE)
  # theme(legend.position = "bottom")
  # browser()
  var_table <- gridExtra::tableGrob(
    pstats %>%
      dplyr::select(component, prop_var) %>%
      dplyr::mutate(prop_var = formatC(prop_var * 100, digits = 3, format = "fg")) %>%
      tidyr::spread(component, prop_var),
    rows = NULL,
    theme = gridExtra::ttheme_default(base_size = 6)
  )
  patchwork::wrap_plots(p_plot, var_table, heights = c(5, 1), ncol = 1)
}

## Save all combinations of PC1-PC5

pca_plot_pc_param <- function(matrix, x, y, col_annotation, facet = NULL) {
  pca_plot(
    matrix,
    col_annotation %>% arrange(condition),
    aes(!!ensym(x), !!ensym(y), fill = as.factor(ERKi), size = as.factor(Time), color = as.character(DOX)),
    center = FALSE, scale = FALSE,
    extra_layers = list(
      geom_point(shape = 21, stroke = 2),
      scale_color_manual(values = c("0" = "#000000", "1" = "#00000000")),
      scale_fill_viridis_d(),
      guides(fill = guide_legend(override.aes = list(color = "#00000000"))),
      # scale_size_manual(values = c("0" = 1, "1" = 2, "2" = 3, "4" = 4, "8" = 5, "16" = 6, "24" = 7)),
      # if (!is.null(facet)) facet_wrap(vars(!!sym(facet))) else NULL,
      facet_grid(rows = vars(Time), cols = vars(ERKi)),
      ggrepel::geom_text_repel(aes(label = replicate), size = 4, color = "black", max.overlaps = Inf)
      # geom_text(aes(label = Repeat), color = "black")
    )
  )
}

x <- pca_plot_pc_param(
  t(counts_vst),
  PC1, PC2,
  meta %>%
    mutate(condition = Sample_ID, across(replicate, as.character)),
  facet = "Time"
)

pca_plots <- combn(paste0("PC", 1:6), 2) %>%
  t() %>%
  `colnames<-`(c("x", "y")) %>%
  as_tibble() %>%
  crossing(
    tibble(facet = list(NULL, "Time", "ERKi"))
  ) %>%
  mutate(
    plot = pmap(
      .,
      pca_plot_pc_param,
      col_annotation = col_annotation
    )
  )

# Cross-correlation -----------------------------------------------------------
###############################################################################T

correlations <- deseq_pairwise %>%
  counts(normalized = TRUE) %>%
  cor() %>%
  as.table() %>%
  as.data.frame() %>%
  as_tibble() %>%
  magrittr::set_colnames(c("Sample_1", "Sample_2", "correlation"))

meta_correlation <- meta %>%
  arrange(
    ERKi, Time, DOX, Repeat, replicate
  ) %>%
  group_by(condition, Repeat) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    Sample_ID,
    sample_name = paste0(condition, "_", replicate) %>%
      fct_inorder(),
    condition
  )

correlation_heatmap_data <- correlations %>%
  inner_join(
    meta_correlation %>%
      dplyr::rename(sample_name_1 = sample_name, condition_1 = condition),
    by = c("Sample_1" = "Sample_ID")
  ) %>%
  inner_join(
    meta_correlation %>%
      dplyr::rename(sample_name_2 = sample_name, condition_2 = condition),
    by = c("Sample_2" = "Sample_ID")
  )

correlation_heatmap <- correlation_heatmap_data %>%
  ggplot(
    aes(sample_name_1, sample_name_2, fill = correlation)
  ) +
  geom_tile() +
  scale_fill_viridis_c(direction = -1, limits = c(0.90, 1), oob = scales::squish) +
  # scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(0.90, 1), oob = scales::squish) +
  geom_rect(
    aes(
      xmin = xmin, xmax = xmax,
      ymin = xmin, ymax = xmax
    ),
    inherit.aes = FALSE,
    color = "black",
    fill = NA,
    data = correlation_heatmap_data %>%
      mutate(across(starts_with("sample_name"), as.integer)) %>%
      group_by(condition_1, condition_2) %>%
      summarize(
        xmin = head(sample_name_1, n = 1) - 0.5, xmax = tail(sample_name_2, n = 1) + 0.5,
        ymin = head(sample_name_2, n = 1) - 0.5, ymax = tail(sample_name_1, n = 1) + 0.5,
      ) %>%
      ungroup()
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("correlation_heatmap.pdf", correlation_heatmap, width = 15, height = 12)

correlations <- deseq_pairwise %>%
  counts(normalized = TRUE) %>%
  cor() %>%
  as.table() %>%
  as.data.frame() %>%
  as_tibble() %>%
  magrittr::set_colnames(c("Sample_1", "Sample_2", "correlation"))

meta_correlation <- meta %>%
  group_by(condition, Repeat) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(
    ERKi, Time, DOX, Repeat, replicate
  ) %>%
  transmute(
    condition,
    Sample_ID,
    sample_name = paste0(condition, "_", Repeat) %>%
      fct_inorder(),
    Repeat = paste0("Repeat_", Repeat)
  ) %>% {
    # browser()
    df <- .
    pivot_wider(
      .,
      id_cols = "condition",
      names_from = "Repeat",
      values_from = c("sample_name", "Sample_ID")
    ) %>%
      tidyr::expand(
        Sample_ID_Repeat_1, Sample_ID_Repeat_2
      ) %>%
      inner_join(
        rename_with(df, ~paste0(.x, "_Repeat_1")),
        by = "Sample_ID_Repeat_1"
      ) %>%
      inner_join(
        rename_with(df, ~paste0(.x, "_Repeat_2")),
        by = c("Sample_ID_Repeat_2")
      )
  }



correlation_heatmap_data <- correlations %>%
  inner_join(
    meta_correlation,
    by = c("Sample_1" = "Sample_ID_Repeat_1", "Sample_2" = "Sample_ID_Repeat_2")
  )

correlation_heatmap <- correlation_heatmap_data %>%
  ggplot(
    aes(sample_name_Repeat_1, sample_name_Repeat_2, fill = correlation)
  ) +
  geom_tile() +
  scale_fill_viridis_c(direction = -1, limits = c(0.90, 1), oob = scales::squish) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("correlation_heatmap_pairwise.pdf", correlation_heatmap, width = 12, height = 10)