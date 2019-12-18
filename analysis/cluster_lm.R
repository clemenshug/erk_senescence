library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(polynom)
library(broom)
library(furrr)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("clustering", "linear_model_clustering")
dir.create(wd, recursive = TRUE, showWarnings = FALSE)

# set directories, import files ------------------------------------------------
###############################################################################T

pairwise_lfc <- syn("syn21432183") %>%
  read_csv()
meta <- syn("syn21432975") %>%
  read_csv()
pca_values <- syn("syn21444456") %>%
  read_csv()
surface_fit <- syn("syn21444486") %>%
  read_csv()

condition_meta <- meta %>%
  distinct(condition, ERKi, Time, DOX)

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

# Find maximum of PC1-2 curve --------------------------------------------------
###############################################################################T
#
poly_fit <- lm(PC2 ~ poly(PC1, 2, raw = TRUE), data = pca_values)

poly_fun <- polynomial(coef = coef(poly_fit))
mid_point <- solve(deriv(poly_fun))

ggplot(pca_values, aes(PC1, PC2)) +
  geom_line() +
  geom_point() +
  geom_smooth(formula = y ~ poly(x, 2), method = "lm") +
  geom_vline(xintercept = mid_point)

# Divide data in high and low erk portion --------------------------------------
###############################################################################T

erk_range_meta <- condition_meta %>%
  filter(condition != "time0dox0conc0", Time == 24) %>%
  left_join(select(pca_values, condition, PC1, PC2), by = "condition") %>%
  mutate(
    erk_range = cut(PC1, c(-Inf, mid_point, Inf), labels = c("low", "high"))
  )

pc_lfc_long <- lfc_long %>%
  # filter(gene_id %in% surface_fit$gene_id) %>%
  inner_join(select(erk_range_meta, condition, PC1, PC2, erk_range), by = "condition")

erk_range_lfc <- pc_lfc_long  %>%
  group_nest(erk_range, gene_id, gene_name)

# Fit linear model to each part ------------------------------------------------
###############################################################################T

p_cutoff <- 0.05
estimate_cutoff <- 0.01

plan(multisession(workers = 4))
erk_range_lfc_fit <- erk_range_lfc %>%
  # crossing(PC_cor = c("PC1", "PC2")) %>%
  # broom::tidy(conf.int = TRUE) %>%
  # filter(term != "(Intercept)") %>%
  # mutate(n = nrow(.x)),
  # mutate(
  #   data = map2(
  #     data, PC_cor,
  #     ~cor.test(.x[[.y]], .x$log2FoldChange, method = "kendall") %>%
  #       tidy()
  #   )
  # ) %>%
  # mutate(
  #   padj = p.adjust(p.value, method = "fdr"),
  #   signed_effect = -log10(padj) * sign(estimate),
  #   class = case_when(
  #     conf.high >= 0.05 ~ "+",
  #     conf.low > -0.05 & conf.high < 0.05 ~ "0",
  #     conf.low <= -0.05 ~ "-",
  #     TRUE ~ "NA"
  #   )
  # )
  mutate(
    linear_model = future_map(
      data,
      ~lm(log2FoldChange ~ PC1, data = .x),
      .progress = TRUE
    ),
    linear_df = map(
      linear_model,
      ~broom::tidy(.x, conf.int = TRUE) %>%
      filter(term != "(Intercept)")
    )
  ) %>%
  unnest(linear_df) %>%
  mutate(
    padj = p.adjust(p.value, method = "fdr"),
    signed_effect = -log10(padj) * sign(estimate),
    # class = cut(signed_effect, c(-Inf, log10(p_cutoff), -log10(p_cutoff), +Inf), labels = c("-", "0", "+"))
    class = case_when(
      signed_effect <= log10(p_cutoff) & estimate < -estimate_cutoff ~ "-",
      signed_effect >= -log10(p_cutoff) & estimate > estimate_cutoff ~ "+",
      TRUE ~ "0"
    )
  )

# Divide genes into classes according to their dynamics ------------------------
###############################################################################T

classes <- tribble(
  ~class, ~direction, ~low, ~high,
  "low_erk", "-", "-", "0",
  "high_erk", "-", "0", "-",
  "bell", "-", "-", "+",
  "full_range", "-", "-", "-",
  "low_erk", "+", "+", "0",
  "high_erk", "+", "0", "+",
  "bell", "+", "+", "-",
  "full_range", "+", "+", "+",
  "no_response", "0", "0", "0",
) %>%
  mutate(class_combined = paste(class, direction, sep = "_"))

erk_range_lfc_classes <- erk_range_lfc_fit %>%
  select(gene_id, gene_name, erk_range, class) %>%
  spread(erk_range, class) %>%
  full_join(classes, by = c("low", "high"))

# Plot representative members of each class ------------------------------------
###############################################################################T

pc_vs_lfc_plot_grid <- function(pca_lfc, genes, aesthetics = aes(PC_value, log2FoldChange, color = PC1), extra_layers = NULL) {
  plot_df <- if (!missing(genes)) {
    pca_lfc %>%
      filter(gene_id %in% genes)
  } else {
    pca_lfc
  }
  plot_df %>%
    ggplot(aesthetics) +
    # geom_point() +
    # facet_grid(vars(gene_name), vars(PC), scales = "free") +
    # geom_point(shape = 21, stroke = 2) +
    geom_path() +
    # scale_color_manual(values = c("0" = "#000000", "1" = "#00000000")) +
    scale_fill_viridis_d() +
    guides(fill = guide_legend(override.aes = list(color = "#00000000"))) +
    scale_size_manual(values = c("1" = 1, "2" = 2, "4" = 3, "8" = 4, "16" = 5, "24" = 6)) +
    guides(color = FALSE) +
    extra_layers
}

plot_cluster_trajectories <- function(
  data, x, y, condition_id, order_id, trajectory_id, facet_x = NULL, facet_y = NULL,
  all_traces = FALSE
) {
  x_quo <- enquo(x)
  y_quo <- enquo(y)
  condition_id_quo <- enquo(condition_id)
  trajectory_id_quo <- enquo(trajectory_id)
  order_id_quo <- enquo(order_id)
  facet_x_quo <- enquo(facet_x)
  facet_y_quo <- enquo(facet_y)
  averages <- data %>%
    group_by(!!order_id_quo, !!condition_id_quo, !!facet_x_quo, !!facet_y_quo) %>%
    summarize_at(vars(!!x_quo, !!y_quo), mean) %>%
    ungroup() %>%
    arrange(!!order_id_quo)
  p <- ggplot(averages, aes(!!x_quo, !!y_quo, color = !!order_id_quo))
  if (all_traces)
    p <- p +
    geom_path(
      aes(group = !!trajectory_id_quo),
      data = arrange(data, !!order_id_quo), size = .8, alpha = 0.1, color = "#CCCCCC"
    )
  p <- p +
    geom_path() +
    geom_point() +
    theme_minimal() +
    scale_color_viridis_c() +
    theme(legend.position = "bottom")
  if (!missing(facet_x) || !missing(facet_y))
    p <- p +
    facet_grid(vars(!!facet_y_quo), vars(!!facet_x_quo), scales = "free")
  p
}

set.seed(42)
lm_cluster_example_traces <- pc_vs_lfc_plot_grid(
  pc_lfc_long %>%
    arrange(PC1) %>%
    mutate(PC_value = PC1) %>%
    inner_join(erk_range_lfc_classes, by = c("gene_id", "gene_name")) %>%
    mutate(wrapping = paste(class, direction, gene_name, sep = "_")),
  erk_range_lfc_classes %>%
    filter(gene_id %in% surface_fit$gene_id) %>%
    group_by(class, direction) %>%
    group_modify(~sample_n(.x, size = 5)) %>%
    pull(gene_id),
  extra_layers = list(theme_minimal(), scale_color_viridis_c(), facet_wrap(vars(wrapping), scales = "free_y"))
)
ggsave(
  file.path(wd, "lm_clustering_example_traces.pdf"),
  lm_cluster_example_traces, height = 8, width = 12
)


lm_cluster_averages <- plot_cluster_trajectories(
  pc_lfc_long %>%
    group_by(gene_id) %>%
    mutate_at(vars(log2FoldChange), scale, scale = TRUE, center = TRUE) %>%
    ungroup() %>%
    inner_join(erk_range_lfc_classes, by = c("gene_id", "gene_name")) %>%
    mutate(PC1_order = PC1) %>%
    gather("PC", "PC_value", PC1, PC2),
  PC_value, log2FoldChange, condition, PC1_order, gene_id, PC, class_combined,
  all_traces = TRUE
)
ggsave(
  file.path(wd, "lm_clustering_average_traces.pdf"),
  lm_cluster_averages, height = 10, width = 3
)

write_csv(
  erk_range_lfc_classes,
  file.path(wd, "lm_clustering_classes.csv")
)

lm_estimates_scatter_plot <- erk_range_lfc_fit %>%
  select(erk_range,  gene_id, signed_effect) %>%
  spread(erk_range, signed_effect) %>%
  inner_join(select(erk_range_lfc_classes, -low, -high), by = "gene_id") %>%
  # filter(class != "no_response") %>%
  filter(gene_id %in% surface_fit$gene_id) %>%
  ggplot(aes(low, high, color = class, shape = direction)) +
    geom_point() +
    # facet_wrap(vars(PC_cor)) +
    geom_density_2d()
    # scale_x_log10() +
    # scale_y_log10() +
    # geom_smooth(method = "lm", formula = I(y - intercept) ~ 0 + x)
ggsave(
  file.path(wd, "lm_estimates_scatter.pdf"),
  lm_estimates_scatter_plot, height = 5, width = 6
)


# Fit different models to the entire range of ERK at once ----------------------
###############################################################################T

erk_range_lfc_model_fit <- pc_lfc_long %>%
  group_nest(gene_id, gene_name) %>%
  mutate(
    linear_model = future_map(
      data,
      ~lm(log2FoldChange ~ PC1, data = .x),
      .progress = TRUE
    ),
    quad_model = future_map(
      data,
      ~lm(log2FoldChange ~ poly(PC1, 2), data = .x),
      .progress = TRUE
    ),
    exp_model = future_map(
      data,
      ~possibly(nls, NULL)(log2FoldChange ~ SSfpl(PC1, A, B, xmid, scal), data = .x),
      .progress = TRUE
    )
  )

erk_range_lfc_model_fit_aic <- erk_range_lfc_model_fit %>%
  mutate_at(vars(ends_with("_model")), map_dbl, possibly(AIC, NA_real_), k = log(12))

erk_range_lfc_model_fit_aic_best_fit <- erk_range_lfc_model_fit_aic %>%
  select(-data) %>%
  gather("model", "aic", -starts_with("gene_")) %>%
  group_by(gene_id, gene_name) %>%
  arrange(aic) %>%
  summarize(model = head(model, 1)) %>%
  ungroup()


set.seed(42)
pc_vs_lfc_plot_grid(
  pc_lfc_long %>%
    arrange(PC1) %>%
    inner_join(erk_range_lfc_model_fit_aic_best_fit, by = c("gene_id", "gene_name")) %>%
    mutate(PC_value = PC1) %>%
    # inner_join(erk_range_lfc_classes, by = c("gene_id", "gene_name")) %>%
    mutate(wrapping = paste(model, gene_name, sep = "_")),
  erk_range_lfc_model_fit_aic_best_fit %>%
    group_by(model) %>%
    group_modify(~sample_n(.x, 10)) %>%
    pull(gene_id),
  extra_layers = list(theme_minimal(), scale_color_viridis_c(), facet_wrap(vars(wrapping), scales = "free_y"))
)


# Do thresholding on moderated slope estimates using ashr ----------------------
###############################################################################T

adjusted_effect_cutoff <- 0.01

erk_range_lfc_fit_ashr <- erk_range_lfc_fit %>%
  bind_cols(
    with(
      .,
      ashr::ash(estimate, std.error, df = 4, mixcompdist = "uniform", method = "fdr")$result
    )
  ) %>%
  mutate(
    class = cut(PosteriorMean, c(-Inf, -adjusted_effect_cutoff, adjusted_effect_cutoff, +Inf), labels = c("-", "0", "+"))
  )

# Clustering genes based on slopes for low and high ERK ------------------------
###############################################################################T

library(mclust)

clPairs(select(erk_range_lfc_fit_wide, -starts_with("gene")))

mclust_bic <- mclustBIC(select(erk_range_lfc_fit_wide, -starts_with("gene")))


mclust_fit <-  Mclust(
  select(erk_range_lfc_fit_wide, -starts_with("gene")),
  G = 8, modelNames = "VVE"
)

plot(mclust_fit, what = "classification")
