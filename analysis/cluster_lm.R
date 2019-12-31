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

syn_lm_clusters <- "syn21448567"

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

estimate_asym_prior <- function(data) {
  l <- lm(y ~ x, data = data)
  coefs <- coef(l)
  data <- mutate(data, diffs = residuals(l))
  diffs_poly <- lm(diffs ~ poly(x, 2), data = data)
  coefs_poly <- coef(diffs_poly)
  asym <- if_else(coefs_poly[[3]] > 0, min(data[["y"]]), max(data[["y"]]))
  l_abs_log <- lm(labs_y ~ x, data = mutate(data, labs_y = log(abs(y - asym) + 1)))
  # browser()
  c(
    a = coef(l_abs_log)[[2]],
    b = asym,
    c = if_else(coefs_poly[[3]] > 0, min(data[["y"]]) + .1, max(data[["y"]]) - 0.1)
  )
}


erk_range_lfc_model_fit <- pc_lfc_long %>%
  left_join(condition_meta, by = "condition") %>%
  filter(Time == 24) %>%
  group_nest(gene_id, gene_name) %>%
  mutate(
    # linear_model = map(
    #   data,
    #   ~lm(log2FoldChange ~ PC1, data = .x)
    #   # .progress = TRUE
    # ),
    # quad_model =  map(
    #   data,
    #   ~lm(log2FoldChange ~ poly(PC1, 2), data = .x)
    #   # .progress = TRUE
    # ),
    # exp_model = future_map(
    #   data,
    #   ~possibly(nls, NULL)(log2FoldChange ~ SSfpl(PC1, A, B, xmid, scal), data = .x),
    #   .progress = TRUE
    # )
    linear_model = map(
      data,
      ~possibly(nls, NULL)(log2FoldChange ~ a*PC1 + x0, data = .x, start = c(a = 1, x0 = 0))
      # .progress = TRUE
    ),
    quad_model =  map(
      data,
      ~possibly(nls, NULL)(log2FoldChange ~ a * PC1**2 + b * PC1 + c, data = .x, start = c(a = 1, b = 1, c = 0))
      # .progress = TRUE
    ),
    # asymptotic_model = map(
    #   data,
    #   function(d) {
    #     # browser()
    #     possibly(nls, NULL)(log2FoldChange ~ aomisc::NLS.asymReg(PC1, init, m, plateau), data = d)
    #   }
    #   # .progress = TRUE
    # )
    # asymptotic_model = map(
    #   data,
    #   function(d) {
    #     # browser()
    #     possibly(nls, NULL)(log2FoldChange ~ b - ((b - c) * exp(-a * PC1)), data = d, start = estimate_asym_prior(tibble(x = d$PC1, y = d$log2FoldChange)))
    #   }
    #   # .progress = TRUE
    # )
    # asymptotic_model2 = map(
    #   data,
    #   function(d) {
    #     library(aomisc)
    #     library(drc)
    #     drm(log2FoldChange ~ PC1, fct = DRC.asymReg(), data = d)
    #   }
    #   # .progress = TRUE
    # )
  ) %>%
  gather("model", "model_object", ends_with("_model")) %>%
  mutate(
    aic = map_dbl(model_object, possibly(AIC, NA_real_), k = log(12)),
    model_df = map(model_object, possibly(broom::tidy,  NULL)),
    coefs = map(model_object, coef),
    # p = map_dbl(
    #   coefs,
    #   ~if_else(
    #     is.null(.x),
    #     NA_real_,
    #     psych::harmonic.mean(filter(.x, !term %in% c("c", "b", "x0", "R0", "Asym"))$p.value)
    #   )
    # ),
    p = map2_dbl(
      model_df, model,
      function(co, m) {
        if (is.null(co) || nrow(co) == 0)
          return(NA_real_)
        co %>%
          filter(term %in% c("a")) %>%
          chuck("p.value", 1)
      }
    ),
    padj = p.adjust(p, method = "BH")
  ) %>%
  mutate(
    data = pmap(
      list(data, model, model_object),
      function(d, m, mo) {
        if (is.null(mo))
          return(mutate(d, yfit = NA_real_))
        c <- coef(mo)
        switch(
          m,
          linear_model = d %>%
            mutate(yfit = c[["a"]] * PC1 + c[["x0"]]),
          quad_model = d %>%
            mutate(yfit = c[["a"]] * PC1**2 + c[["b"]] * PC1 + c[["c"]]),
          asymptotic_model = d %>%
            mutate(yfit = c[["b"]] - ((c[["b"]] - c[["c"]]) * exp(-c[["a"]] * PC1)))
        )
      }
    )
  )


erk_range_lfc_model_fit_plotting <- erk_range_lfc_model_fit %>%
  mutate(
    coef_str = map_chr(
      model_df,
      function(d) {
        if (is.null(d) || !is.numeric(d$estimate))
          return("")
        paste(
          d$term,
          formatC(signif(d$estimate, digits = 2), digits = 2, format="fg", flag="#"),
          formatC(d$p.value, digits = 2, format="e"),
          sep = " ", collapse = "\n"
        )
      }
    )
  )


erk_range_lfc_model_best_fit_aic <- erk_range_lfc_model_fit %>%
  drop_na(aic) %>%
  group_by(gene_id, gene_name) %>%
  arrange(aic, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

erk_range_lfc_model_best_fit_p <- erk_range_lfc_model_fit %>%
  drop_na(p) %>%
  filter(padj < 0.05) %>%
  group_by(gene_id, gene_name) %>%
  arrange(p, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

erk_range_lfc_model_best_fit_aic %>%
  count(model)

set.seed(42)
pc_vs_lfc_plot_grid(
  erk_range_lfc_model_fit_plotting %>%
    dplyr::select(gene_name, gene_id, data, model, coef_str) %>%
    semi_join(erk_range_lfc_model_best_fit_aic, by = c("gene_id", "gene_name", "model")) %>%
    unnest(data) %>%
    arrange(PC1) %>%
    filter(Time == 24) %>%
    mutate(PC_value = PC1) %>%
    # inner_join(erk_range_lfc_classes, by = c("gene_id", "gene_name")) %>%
    mutate(wrapping = paste(model, gene_name, sep = "_")),
  erk_range_lfc_model_best_fit_aic %>%
    group_by(model) %>%
    group_modify(~sample_n(.x, 10, replace = TRUE)) %>%
    pull(gene_id),
  # erk_range_lfc_model_best_fit_p %>%
  #   filter(padj > 0.05) %>%
  #   pull(gene_id),
  extra_layers = list(
    theme_minimal(),
    scale_color_viridis_c(),
    facet_wrap(vars(wrapping), scales = "free_y"),
    geom_line(aes(y = yfit), color = "black"),
    geom_text(aes(label =  coef_str), color = "black", x = -Inf, y = Inf, hjust = 0, vjust = 1)
  )
)

classify_full_range <- function(model, coefs) {
  if (model == "linear_model") {
    class <- "full_range"
    direction <- if_else(coefs[["a"]] > 0, "+", "-")
  } else if (model == "quad_model") {
    zc <- -coefs[["b"]] / (2 * coefs[["a"]])
    direction <- if_else(coefs[["a"]] > 0, "-", "+")
    if (abs(zc) < 25) {
      class <- "bell"

    }
    else {
      class <- if_else(zc > 25, "low_erk", "high_erk")
    }
  } else if (model == "asymptotic_model") {

  }
}

full_range_clustering_classes <- erk_range_lfc_model_best_fit_p %>%
  mutate(
    class = case_when(
      model == "linear_model" & map_chr(coefs, "a") > 0 ~ ""
    )
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

# Store to synapse -------------------------------------------------------------
###############################################################################T


lm_clustering_activity <- Activity(
  "Clustering ERK response genes based on their expression profile with increasing ERK concentration",
  used = c(
    "syn21432183",
    "syn21432975",
    "syn21444456",
    "syn21444486"
  ),
  executed = "https://github.com/clemenshug/erk_senescence/blob/master/analysis/lm_clusters.R"
)

c(
  file.path(wd, "lm_clustering_classes.csv")
) %>%
  synStoreMany(syn_lm_clusters, activity = lm_clustering_activity)

