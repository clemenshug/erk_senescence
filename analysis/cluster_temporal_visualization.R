library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(pheatmap)
library(RColorBrewer)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

wd <- here("temporal_ordering")
dir.create(wd, showWarnings = FALSE)

# set directories, import files ------------------------------------------------
###############################################################################T

surface_fit <- syn("syn21444486") %>%
  read_csv()

temporal_ordering <- syn("syn25471789") %>%
  read_csv()

temporal_lfc <- syn("syn25471787") %>%
  read_csv(col_types = list(ERKi = col_character()))

temporal_lfc_rescaled <- syn("syn25471791") %>%
  read_csv()

temporal_lfc_all_stats <- syn("syn25471794") %>%
  read_csv(
    col_types = list(ERKi = col_character())
  )

meta <- syn("syn21432975") %>%
  read_csv()

dose_clusters <- syn("syn21576614") %>%
  read_csv()

condition_meta <- meta %>%
  distinct(condition, ERKi, Time, DOX)

pairwise_lfc <- syn("syn21432183") %>%
  read_csv()

pairwise_padj <- syn("syn21432184") %>%
  read_csv()

padj_long <- pairwise_padj  %>%
  mutate_at(vars(starts_with("time")), replace_na, 0) %>%
  gather("condition", "padj", -gene_name, -gene_id)

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

cluster_names <- syn("syn21567677") %>%
  read_csv() %>%
  mutate(
    class_name = factor(class_name, levels = c("full ERK", "low ERK", "high ERK", "bell")),
    direction = factor(direction, levels = c("upregulated", "downregulated"))
  )

theme_set(
  theme_minimal(base_family = "Helvetica") +
    theme(
      strip.text = element_text(size = rel(1)),
      text = element_text(face = "plain")
    )
)

# Visualize interaction of dose and time clusters ------------------------------
###############################################################################T

time_dose_clusters <- temporal_ordering %>%
  filter(directed == "directed", ERKi %in% c("all", "0", "1000")) %>%
  select(ERKi, gene_id, gene_name, mid_induction_time) %>%
  # pivot_longer(ends_with("induction"), names_to = "method", values_to = "time") %>%
  inner_join(
    dose_clusters %>%
      select(
        gene_id, gene_name, dose_cluster = k_medoids
      ),
    by = c("gene_id", "gene_name")
  ) %>%
  inner_join(
    cluster_names %>%
      select(dose_cluster = class_combined, dose_cluster_name_combined = combined,
             dose_cluster_direction = direction, dose_cluster_name = class_name),
    by = c("dose_cluster")
  ) %>%
  arrange(mid_induction_time) %>%
  mutate(across(mid_induction_time, . %>% as.character() %>% fct_inorder()))

temporal_lfc_rescaled_with_clusters <- temporal_lfc_all_stats %>%
  inner_join(
    dose_clusters %>%
      select(
        gene_id, gene_name, dose_cluster = k_medoids
      ),
    by = c("gene_id", "gene_name")
  ) %>%
  inner_join(
    cluster_names %>%
      select(dose_cluster = class_combined, dose_cluster_name_combined = combined,
             dose_cluster_direction = direction, dose_cluster_name = class_name),
    by = c("dose_cluster")
  )

time_dose_clusters_plot <- time_dose_clusters %>%
  ggplot(
    aes(mid_induction_time, fill = dose_cluster_name)
  ) +
    geom_bar(
      position = "fill"
    ) +
    facet_grid(vars(dose_cluster_direction), vars(ERKi))

cowplot::ggsave2(
  file.path(wd, "time_dose_clusters_barplot.pdf"),
  time_dose_clusters_plot
)

time_dose_clusters_plot <- time_dose_clusters %>%
  filter(ERKi == "all") %>%
  ggplot(
    aes(mid_induction_time, fill = dose_cluster_name)
  ) +
  geom_bar(
    position = "fill"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  facet_wrap(~dose_cluster_direction) +
  labs(x = "Induction time", y = "Proportion", fill = "Cluster")

cowplot::ggsave2(
  file.path(wd, "time_dose_clusters_barplot_only_all.pdf"),
  time_dose_clusters_plot,
  width = 5, height = 4
)


time_dose_clusters <- temporal_ordering %>%
  filter(directed == "absolute", ERKi %in% c("all", "0", "1000")) %>%
  select(ERKi, gene_id, max_induction_time, mid_induction_time) %>%
  pivot_longer(ends_with("induction_time"), names_to = "method", values_to = "time") %>%
  inner_join(
    dose_clusters %>%
      select(gene_id, dose_class = k_medoids) %>%
      filter(dose_class != "no_response_0") %>%
      mutate(
        direction = str_extract(dose_class, "[+-]"),
        dose_class = str_replace(dose_class, "_[+-]", "")
      )
  ) %>%
  arrange(time) %>%
  mutate(across(time, . %>% as.character() %>% fct_inorder()))



time_dose_clusters_plot <- time_dose_clusters %>%
  ggplot(
    aes(dose_class, fill = time, label = stat(count))
  ) +
  geom_bar(
    position = "fill",
    # position = "dodge"
  ) +
  scale_fill_viridis_d(direction = -1) +
  # facet_wrap(vars(method))
  # geom_text(
  #   position = position_dodge(),
  #   stat = "count"
  # ) +
  facet_grid(vars(method), vars(ERKi), scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

cowplot::ggsave2(
  file.path(wd, "time_dose_clusters_barplot_by_class_no_direction.pdf"),
  time_dose_clusters_plot, width = 8, height = 7
)


time_dose_clusters <- temporal_ordering %>%
  filter(directed == "directed", ERKi %in% c("all", "0", "1000")) %>%
  select(ERKi, gene_id, max_induction_time, mid_induction_time) %>%
  pivot_longer(ends_with("induction_time"), names_to = "method", values_to = "time") %>%
  inner_join(
    dose_clusters %>%
      select(gene_id, dose_class = k_medoids) %>%
      filter(dose_class != "no_response_0")
  ) %>%
  semi_join(
    padj_long %>%
      filter(padj < 0.05)
  ) %>%
  arrange(time) %>%
  mutate(across(time, . %>% as.character() %>% fct_inorder()))

time_dose_clusters_plot <- time_dose_clusters %>%
  filter(method == "mid_induction_time") %>%
  inner_join(
    cluster_names,
    by = c("dose_class" = "class_combined")
  ) %>%
  ggplot(
    aes(class_name, fill = time, label = stat(count))
  ) +
  geom_bar(
    position = "fill",
    # position = "dodge"
  ) +
  scale_fill_viridis_d(direction = -1) +
  # facet_wrap(vars(method))
  # geom_text(
  #   position = position_dodge(),
  #   stat = "count"
  # ) +
  facet_grid(vars(direction), vars(ERKi)) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica"),
    # axis.title = element_text(face = "bold"),
    # axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)
  ) +
  labs(x = "Cluster", y = "Proportion")

cowplot::ggsave2(
  file.path(wd, "time_dose_clusters_barplot_by_class_with_direction.pdf"),
  time_dose_clusters_plot, width = 8, height = 7
)


time_dose_clusters_plot <- time_dose_clusters %>%
  filter(method == "mid_induction_time", ERKi == "all") %>%
  inner_join(
    cluster_names,
    by = c("dose_class" = "class_combined")
  ) %>%
  ggplot(
    aes(class_name, fill = time, label = stat(count))
  ) +
  geom_bar(
    position = "fill",
    # position = "dodge"
  ) +
  scale_fill_viridis_d(direction = -1) +
  # facet_wrap(vars(method))
  # geom_text(
  #   position = position_dodge(),
  #   stat = "count"
  # ) +
  facet_wrap(~direction, nrow = 1L) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica"),
    # axis.title = element_text(face = "bold"),
    # axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)
  ) +
  labs(x = "Dose-response behavior", y = "Proportion", fill = "Time")

cowplot::ggsave2(
  file.path(wd, "time_dose_clusters_barplot_by_class_with_direction_only_all_.pdf"),
  time_dose_clusters_plot, width = 4, height = 3.5
)


time_dose_clusters_plot <- temporal_lfc_rescaled_with_clusters %>%
  mutate(
    across(
      c(mid_induction_time),
      . %>%
        factor(., levels = as.character(sort(unique(.))))
    )
  ) %>%
  filter(ERKi == "all", directed == "absolute", log2FoldChange_range_erki > 1) %>%
  ggplot(
    aes(dose_cluster_name, fill = mid_induction_time, label = stat(count))
  ) +
  geom_bar(
    position = "fill",
    # position = "dodge"
  ) +
  scale_fill_viridis_d(direction = -1) +
  # facet_wrap(vars(method))
  # geom_text(
  #   position = position_dodge(),
  #   stat = "count"
  # ) +
  facet_wrap(~dose_cluster_direction, nrow = 1L) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica"),
    # axis.title = element_text(face = "bold"),
    # axis.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)
  ) +
  labs(x = "Cluster", y = "Proportion", fill = "Time")

time_dose_clusters_plot <- temporal_lfc_rescaled_with_clusters %>%
  mutate(
    across(
      c(mid_induction_time),
      . %>%
        factor(., levels = as.character(sort(unique(.))))
    )
  ) %>%
  filter(ERKi == "all", directed == "absolute", log2FoldChange_range_erki > 1) %>%
  ggplot(
    aes(mid_induction_time, fill = dose_cluster_name)
  ) +
  geom_bar(
    position = "fill"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  facet_wrap(~dose_cluster_direction) +
  labs(x = "Induction time", y = "Proportion", fill = "Dose-response\nbehavior")

cowplot::ggsave2(
  file.path(wd, "time_dose_clusters_barplot_only_all_only_sig.pdf"),
  time_dose_clusters_plot,
  width = 5, height = 4
)

time_dose_clusters_plot <- temporal_lfc_rescaled_with_clusters %>%
  mutate(
    across(
      c(mid_induction_time),
      . %>%
        factor(., levels = as.character(sort(unique(.))))
    )
  ) %>%
  filter(ERKi == "all", directed == "absolute", log2FoldChange_range_erki > 1) %>%
  ggplot(
    aes(mid_induction_time, fill = dose_cluster_name)
  ) +
  geom_bar(
    position = "fill"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme_minimal() +
  # facet_wrap(~dose_cluster_direction) +
  labs(x = "Induction time", y = "Proportion", fill = "Cluster")

time_dose_clusters_plot <- temporal_lfc_rescaled_with_clusters %>%
  filter(ERKi == "all", directed == "absolute") %>%
  distinct(gene_id, mid_induction_time, dose_cluster_name, dose_cluster_direction) %>%
  mutate(
    mid_induction_time = factor(
      mid_induction_time,
      levels = as.character(sort(unique(mid_induction_time)))
    )
  ) %>%
  ggplot(aes(
    x = mid_induction_time, fill = dose_cluster_name
  )) +
    geom_bar(aes(y = stat(prop), group = dose_cluster_name), position = "dodge") +
    labs(x = "Induction time", y = "Proportion", fill = "Dose-response\nbehavior") +
    scale_fill_brewer(palette = "Set2")
    # facet_wrap(~dose_cluster_direction)

cowplot::ggsave2(
  file.path(wd, "time_dose_clusters_barplot_by_time_only_all_proportion_per_cluster.pdf"),
  time_dose_clusters_plot,
  width = 5, height = 3
)

time_dose_clusters_plot <- temporal_lfc_rescaled_with_clusters %>%
  filter(ERKi == "all", directed == "absolute") %>%
  distinct(gene_id, mid_induction_time, dose_cluster_name, dose_cluster_direction) %>%
  mutate(
    mid_induction_time = factor(
      mid_induction_time,
      levels = as.character(sort(unique(mid_induction_time)))
    )
  ) %>%
  ggplot(aes(
    x = mid_induction_time, fill = dose_cluster_name
  )) +
  geom_bar(aes(group = dose_cluster_name), position = "dodge") +
  labs(x = "Induction time", y = "Count", fill = "Dose-response\nbehavior") +
  scale_fill_brewer(palette = "Set2")
# facet_wrap(~dose_cluster_direction)

# cowplot::ggsave2(
ggsave(
  file.path(wd, "time_dose_clusters_barplot_by_time_only_all_absolute_counts.pdf"),
  time_dose_clusters_plot,
  width = 5, height = 3,
  device = cairo_pdf
)


time_dose_clusters_plot <- temporal_lfc_rescaled_with_clusters %>%
  filter(ERKi == "all", directed == "absolute") %>%
  distinct(gene_id, mid_induction_time, dose_cluster_name, dose_cluster_direction) %>%
  mutate(
    mid_induction_time = factor(
      mid_induction_time,
      levels = as.character(sort(unique(mid_induction_time)))
    )
  ) %>%
  ggplot(aes(
    x = mid_induction_time, fill = dose_cluster_name
  )) +
  geom_bar(aes(y = stat(prop), group = dose_cluster_name), position = "dodge") +
  labs(x = "Induction time", y = "Proportion", fill = "Dose-response\nbehavior") +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~dose_cluster_direction, ncol = 1)

cowplot::ggsave2(
  file.path(wd, "time_dose_clusters_barplot_by_time_only_all_proportion_per_cluster_separated.pdf"),
  time_dose_clusters_plot,
  width = 5, height = 3
)


time_dose_cluster_independence_test <- temporal_lfc_rescaled_with_clusters %>%
  mutate(
    Time_rank = factor(Time, levels = as.character(sort(unique(Time)))) %>%
      as.numeric()
  ) %>%
  group_nest(directed, ERKi) %>%
  mutate(
    data = map(
      data,
      function(df) {
        coin::kruskal_test(
          mid_induction_time ~ dose_cluster_name,
          df %>%
            distinct(gene_id, mid_induction_time, dose_cluster_name),
          conf.int = TRUE,
          ties.method	= "average-scores"
        )
      }
    )
  )

time_dose_cluster_independence_test <- temporal_lfc_rescaled_with_clusters %>%
  mutate(
    Time_rank = factor(Time, levels = as.character(sort(unique(Time)))) %>%
      as.numeric()
  ) %>%
  distinct(directed, ERKi, dose_cluster_name, mid_induction_time, gene_id) %>%
  group_nest(directed, ERKi) %>%
  mutate(
    data = pmap(
      .,
      function(dose_cluster_name, data, ...) {
        browser()
        data %>%
          bind_rows(
            . %>%
              filter
          )
      }
    )
  )

library(lmPerm)

time_dose_cluster_independence_test <- temporal_lfc_rescaled_with_clusters %>%
  mutate(
    mid_induction_time_rank = factor(mid_induction_time, levels = as.character(sort(unique(mid_induction_time)))) %>%
      as.numeric()
  ) %>%
  distinct(directed, ERKi, dose_cluster_name, mid_induction_time_rank, mid_induction_time, gene_id) %>%
  group_nest(directed, ERKi) %>%
  filter(directed == "absolute", ERKi == "all") %>%
  mutate(
    data = pmap(
      .,
      function(dose_cluster_name, data, ...) {
        # browser()
        # model <- nlme::lme(
        #   # mid_induction_time ~ dose_cluster_name,
        #   mid_induction_time ~ 1,
        #   random = -1 | dose_cluster_name,
        #   data = data
        # )
        # PT <- rcompanion::pairwisePermutationTest(
        #   mid_induction_time ~ dose_cluster_name,
        #   data = data, method = "fdr"
        # )
        # model <- data %>%
        #   specify(mid_induction_time ~ dose_cluster_name) %>%
        #   hypothesise(null = "independence") %>%
        #   generate(reps = 1000, type = "permute")
        # browser()
        # model <- lmPerm::aovp(mid_induction_time ~ dose_cluster_name, data = data)
        bs <- data %>%
          specify(mid_induction_time ~ dose_cluster_name) %>%
          hypothesise(null = "independence") %>%
          generate(reps = 1000, type = "bootstrap")
        bs_stats <- bs %>%
          group_by(dose_cluster_name, replicate) %>%
          summarize(m = mean(mid_induction_time), .groups = "drop") %>%
          group_by(dose_cluster_name) %>%
          summarize(mm = mean(m), ci_low = quantile(m, 0.05), ci_high = quantile(m, 0.95), .group = "drop")
        dunn <- dunn.test::dunn.test(
          x = data$mid_induction_time,
          g = data$dose_cluster_name,
          method = "bh"
        )
        model <- aov(mid_induction_time ~ dose_cluster_name, data = data)
        tukey = TukeyHSD(model)
        list(
          bs = bs,
          bs_stats = bs_stats,
          model = model,
          dunn = dunn,
          tukey = tukey
        )
      }
    )
  )

time_dose_cluster_mean_induction_plot <- ggplot(
  time_dose_cluster_independence_test$data[[1]]$bs %>%
    group_by(dose_cluster_name, replicate) %>%
    summarize(m = mean(mid_induction_time), .groups = "drop"),
  aes(dose_cluster_name, m, fill = dose_cluster_name)
) +
  geom_violin(
    color = "black",
    draw_quantiles = c(0.05, 0.5, 0.95),
    show.legend = FALSE
  ) +
  scale_fill_brewer(palette = "Set2") +
  # geom_boxplot() +
  lims(y = c(0, NA)) +
  labs(x = "Dose-response behavior", y = "Bootstrapped mean induction time")

cowplot::ggsave2(
  file.path(wd, "time_dose_clusters_bootstrapped_mean_violin.pdf"),
  time_dose_cluster_mean_induction_plot,
  width = 3, height = 3.5
)

  mutate(
    data = map(
      data,
      function(df) {
        coin::kruskal_test(
          mid_induction_time ~ dose_cluster_name,
          df %>%
            distinct(gene_id, mid_induction_time, dose_cluster_name),
          conf.int = TRUE,
          ties.method	= "average-scores"
        )
      }
    )
  )

# Time trace of example gene all concentrations --------------------------------
###############################################################################T

gene_temporal_lfc <- temporal_lfc_all_stats %>%
    filter(gene_name == "CDKN1B", directed == "absolute") %>% {
      bind_rows(
        filter(., ERKi == "all") %>%
          select(Time, starts_with("log2FoldChange_q")) %>%
          pivot_longer(-Time, names_to = "series", values_to = "log2FoldChange") %>%
          mutate(
            series = if_else(str_detect(series, fixed("q33")), "33rd quantile", "67th quantile")
          ),
        filter(., ERKi != "all") %>%
          select(Time, series = ERKi, log2FoldChange)
      )
    } %>%
    mutate(
      quantile = str_detect(series, fixed("quantile"))
    )

gene_temporal_lfc_plot <- gene_temporal_lfc %>%
  mutate(
    series = factor(
      series, levels = c(unique(meta$ERKi), "33rd quantile", "67th quantile")
    )
  ) %>%
  ggplot(aes(Time, log2FoldChange, color = series)) +
  geom_line(aes(size = quantile, linetype = quantile)) +
  scale_color_manual(
    values = c(
      set_names(viridisLite::viridis(6), c("0", "3.9", "15.6", "62.5", "250", "1000")),
      "33rd quantile" = "blue", "67th quantile" = "red"
    )
  ) +
  scale_size_manual(
    values = c(`TRUE` = 2, `FALSE` = 1)
  ) +
  scale_linetype_manual(
    values = c(`TRUE` = "dashed", `FALSE` = "solid")
  ) +
  guides(size = FALSE, linetype = FALSE) +
  labs(color = "ERKi dose", title = "CDKN1B")

ggsave(
  file.path(wd, "time_series_explanation_CDKN1B.pdf"),
  gene_temporal_lfc_plot,
  width = 5, height = 3
)

x <- temporal_lfc_all_stats %>%
  filter(gene_name == "CDKN2B", directed == "absolute", ERKi != "all", Time == 24)

y <- wtd.quantile(
  x$log2FoldChange, weights = x$log2FoldChange_range_erki, probs = c(0.33, 0.67)
)

y <- Hmisc::wtd.quantile(
  x$log2FoldChange, weights = x$log2FoldChange_range_erki, probs = c(0.33, 0.67)
)

# Count number of genes flipping signs ------------------------------------------
###############################################################################T


n_genes <- temporal_lfc_rescaled_with_clusters %>%
  # filter(ERKi == "0", directed == "directed", log2FoldChange_range_erki > 1) %>%
  mutate(
    pos_neg = sign(log2FoldChange)
  ) %>%
  group_by(gene_id, gene_name) %>%
  summarize(
    pos_neg_sum = sum(pos_neg),
    .groups = "drop"
  )

n_genes_plot <- n_genes %>%
  mutate(
    pos_neg_sum = factor(
      pos_neg_sum,
      levels = as.character(sort(unique(pos_neg_sum)))
    )
  ) %>%
  ggplot(aes(pos_neg_sum)) +
  geom_bar()



# Visualize induction clusters -------------------------------------------------
###############################################################################T

plot_cluster_trajectories <- function(
  data, x, y, order_id, trajectory_id, facet_x = NULL, facet_y = NULL,
  all_traces = FALSE, all_trace_alpha = 0.2
) {
  x_quo <- enquo(x)
  y_quo <- enquo(y)
  order_id_quo <- enquo(order_id)
  trajectory_id_quo <- enquo(trajectory_id)
  facet_x_quo <- enquo(facet_x)
  facet_y_quo <- enquo(facet_y)
  averages <- data %>%
    group_by(!!x_quo, !!facet_x_quo, !!facet_y_quo) %>%
    summarize_at(vars(!!y_quo), mean) %>%
    ungroup() %>%
    arrange(!!order_id_quo)
  p <- ggplot(averages, aes(!!x_quo, !!y_quo, color = !!order_id_quo))
  if (all_traces)
    p <- p +
    geom_path(
      aes(group = !!trajectory_id_quo),
      data = data %>% arrange(!!order_id_quo), size = .8, alpha = all_trace_alpha, color = "#CCCCCC"
    )
  p <- p +
    geom_path() +
    geom_point() +
    # theme_minimal() +
    scale_color_viridis_c() +
    theme(legend.position = "bottom")
  if (!missing(facet_x) && !missing(facet_y)) {
    p <- p +
      facet_grid(vars(!!facet_y_quo), vars(!!facet_x_quo), scales = "free")
  }
  else (!missing(facet_x))
    p <- p +
      facet_wrap(vars(!!facet_x_quo), scales = "free")
  p
}

temporal_ordering_traces <-  plot_cluster_trajectories(
  temporal_lfc_rescaled_with_clusters %>%
    filter(ERKi == "all", directed == "directed", log2FoldChange_range_erki > 1) %>%
    mutate(
      log2FoldChange_rescaled = log2FoldChange_rescaled * if_else(
        max_induction < 0, -1, 1
      )
    ) %>%
    # select(gene_id, Time, dose_cluster_name_combined, dose_cluster_direction,
    #        dose_cluster_name, mid_induction_time, log2FoldChange_rescaled) %>%
    mutate_at(vars(mid_induction_time), . %>% as.character %>% fct_reorder(as.numeric(.))),
  x = Time, y = log2FoldChange_rescaled,
  order_id = Time, trajectory_id = gene_id,
  facet_x = mid_induction_time, facet_y = dose_cluster_name, all_traces = TRUE,
  all_trace_alpha = 0.4
) +
  labs(y = "Expression induction") +
  facet_grid(vars(dose_cluster_name), vars(mid_induction_time))

cowplot::ggsave2(
  file.path(wd, "temporal_ordering_traces_by_dose_cluster_directed.pdf"),
  temporal_ordering_traces, width = 10, height = 6
)



temporal_lfc_rescaled_with_clusters %>%
  mutate(across(log2FoldChange, round, digits = 2)) %>%
  filter(directed == "directed", ERKi == "all", mid_induction_time == 24, dose_cluster_name == "full ERK") %>%
  View()



temporal_ordering_traces <-  plot_cluster_trajectories(
  temporal_lfc_rescaled_with_clusters %>%
    filter(ERKi == "all", directed == "absolute") %>%
    mutate(
      log2FoldChange_rescaled = log2FoldChange_rescaled * if_else(
        max_induction < 0, -1, 1
      )
    ) %>%
    # select(gene_id, Time, dose_cluster_name_combined, dose_cluster_direction,
    #        dose_cluster_name, mid_induction_time, log2FoldChange_rescaled) %>%
    mutate_at(vars(mid_induction_time), . %>% as.character %>% fct_reorder(as.numeric(.))),
  x = Time, y = log2FoldChange_rescaled,
  order_id = Time, trajectory_id = gene_id,
  facet_x = mid_induction_time, facet_y = dose_cluster_name, all_traces = TRUE,
  all_trace_alpha = 0.4
) +
  labs(y = "Expression induction") +
  facet_grid(vars(dose_cluster_name), vars(mid_induction_time))

cowplot::ggsave2(
  file.path(wd, "temporal_ordering_traces_by_dose_cluster_absolute.pdf"),
  temporal_ordering_traces, width = 10, height = 6
)


temporal_ordering_traces <-  plot_cluster_trajectories(
  temporal_lfc_rescaled_with_clusters %>%
    filter(ERKi == "all", directed == "absolute") %>%
    mutate(
      log2FoldChange_rescaled = log2FoldChange_rescaled * if_else(
        max_induction < 0, -1, 1
      )
    ) %>%
    # select(gene_id, Time, dose_cluster_name_combined, dose_cluster_direction,
    #        dose_cluster_name, mid_induction_time, log2FoldChange_rescaled) %>%
    mutate_at(vars(mid_induction_time), . %>% as.character %>% fct_reorder(as.numeric(.))),
  x = Time, y = log2FoldChange_rescaled,
  order_id = Time, trajectory_id = gene_id,
  facet_x = mid_induction_time, facet_y = dose_cluster_direction, all_traces = TRUE,
  all_trace_alpha = 0.4
) +
  labs(y = "Expression induction") +
  facet_grid(vars(dose_cluster_direction), vars(mid_induction_time))


temporal_ordering_traces <-  plot_cluster_trajectories(
  temporal_lfc_rescaled_with_clusters %>%
    filter(ERKi == "all", directed == "absolute") %>%
    mutate(
      log2FoldChange_rescaled = log2FoldChange_rescaled * if_else(
        max_induction < 0, -1, 1
      )
    ) %>%
    # select(gene_id, Time, dose_cluster_name_combined, dose_cluster_direction,
    #        dose_cluster_name, mid_induction_time, log2FoldChange_rescaled) %>%
    mutate_at(vars(mid_induction_time), . %>% as.character %>% fct_reorder(as.numeric(.))),
  x = Time, y = log2FoldChange_rescaled,
  order_id = Time, trajectory_id = gene_id,
  facet_x = mid_induction_time, facet_y = dose_cluster_name_combined, all_traces = TRUE,
  all_trace_alpha = 0.4
) +
  labs(y = "Expression induction") +
  facet_grid(vars(dose_cluster_name_combined), vars(mid_induction_time))



temporal_ordering_traces <-  plot_cluster_trajectories(
  temporal_lfc_rescaled_with_clusters %>%
    filter(ERKi %in% c("0", "1000"), directed == "absolute") %>%
    mutate(
      log2FoldChange_rescaled = log2FoldChange_rescaled * if_else(
        max_induction < 0, -1, 1
      )
    ) %>%
    # select(gene_id, Time, dose_cluster_name_combined, dose_cluster_direction,
    #        dose_cluster_name, mid_induction_time, log2FoldChange_rescaled) %>%
    mutate_at(vars(mid_induction_time), . %>% as.character %>% fct_reorder(as.numeric(.))),
  x = Time, y = log2FoldChange_rescaled,
  order_id = Time, trajectory_id = gene_id,
  facet_x = mid_induction_time, facet_y = dose_cluster_name_combined, all_traces = TRUE,
  all_trace_alpha = 0.4
) +
  labs(y = "Expression induction") +
  facet_grid(vars(dose_cluster_name_combined), vars(mid_induction_time))



temporal_ordering_traces <-  plot_cluster_trajectories(
  temporal_lfc_rescaled_with_clusters %>%
    filter(ERKi %in% c("1000"), directed == "directed", log2FoldChange_range_erki > 1) %>%
    mutate(
      log2FoldChange_rescaled = log2FoldChange_rescaled * if_else(
        max_induction < 0, -1, 1
      )
    ) %>%
    # select(gene_id, Time, dose_cluster_name_combined, dose_cluster_direction,
    #        dose_cluster_name, mid_induction_time, log2FoldChange_rescaled) %>%
    mutate_at(vars(mid_induction_time), . %>% as.character %>% fct_reorder(as.numeric(.))),
  x = Time, y = log2FoldChange_rescaled,
  order_id = Time, trajectory_id = gene_id,
  facet_x = mid_induction_time, facet_y = dose_cluster_name, all_traces = TRUE,
  all_trace_alpha = 0.4
) +
  labs(y = "Expression induction") +
  facet_grid(vars(dose_cluster_name), vars(mid_induction_time))

cowplot::ggsave2(
  file.path(wd, "temporal_ordering_traces_directed_ERKi_1000.pdf"),
  temporal_ordering_traces, width = 10, height = 6
)

temporal_ordering_traces <-  plot_cluster_trajectories(
  temporal_lfc_rescaled_with_clusters %>%
    filter(ERKi %in% c("0"), directed == "directed", log2FoldChange_range_erki > 1) %>%
    mutate(
      log2FoldChange_rescaled = log2FoldChange_rescaled * if_else(
        max_induction < 0, -1, 1
      )
    ) %>%
    # select(gene_id, Time, dose_cluster_name_combined, dose_cluster_direction,
    #        dose_cluster_name, mid_induction_time, log2FoldChange_rescaled) %>%
    mutate_at(vars(mid_induction_time), . %>% as.character %>% fct_reorder(as.numeric(.))),
  x = Time, y = log2FoldChange_rescaled,
  order_id = Time, trajectory_id = gene_id,
  facet_x = mid_induction_time, facet_y = dose_cluster_name, all_traces = TRUE,
  all_trace_alpha = 0.4
) +
  labs(y = "Expression induction") +
  facet_grid(vars(dose_cluster_name), vars(mid_induction_time))

cowplot::ggsave2(
  file.path(wd, "temporal_ordering_traces_directed_ERKi_0.pdf"),
  temporal_ordering_traces, width = 10, height = 6
)


temporal_ordering_traces <-  plot_cluster_trajectories(
  temporal_lfc_rescaled_with_clusters %>%
    filter(ERKi %in% c("all"), directed == "absolute", log2FoldChange_range_erki > 1) %>%
    mutate(
      log2FoldChange_rescaled = log2FoldChange_rescaled * if_else(
        max_induction < 0, -1, 1
      )
    ) %>%
    # select(gene_id, Time, dose_cluster_name_combined, dose_cluster_direction,
    #        dose_cluster_name, mid_induction_time, log2FoldChange_rescaled) %>%
    mutate_at(vars(mid_induction_time), . %>% as.character %>% fct_reorder(as.numeric(.))),
  x = Time, y = log2FoldChange_rescaled,
  order_id = Time, trajectory_id = gene_id,
  facet_x = mid_induction_time, facet_y = dose_cluster_name, all_traces = TRUE,
  all_trace_alpha = 0.4
) +
  labs(y = "Expression induction") +
  facet_grid(vars(dose_cluster_name), vars(mid_induction_time))

cowplot::ggsave2(
  file.path(wd, "temporal_ordering_traces_absolute_ERKi_all.pdf"),
  temporal_ordering_traces, width = 10, height = 6
)

temporal_ordering_traces <-  plot_cluster_trajectories(
  temporal_lfc_rescaled_with_clusters %>%
    filter(ERKi %in% c("all"), directed == "absolute", log2FoldChange_range_erki > 1) %>%
    # select(gene_id, Time, dose_cluster_name_combined, dose_cluster_direction,
    #        dose_cluster_name, mid_induction_time, log2FoldChange_rescaled) %>%
    mutate_at(vars(mid_induction_time), . %>% as.character %>% fct_reorder(as.numeric(.))),
  x = Time, y = log2FoldChange_rescaled,
  order_id = Time, trajectory_id = gene_id,
  facet_x = mid_induction_time, facet_y = dose_cluster_name, all_traces = TRUE,
  all_trace_alpha = 0.4
) +
  labs(y = "Expression induction") +
  facet_grid(vars(dose_cluster_name), vars(mid_induction_time))

cowplot::ggsave2(
  file.path(wd, "temporal_ordering_traces_absolute_ERKi_all_no_flipping.pdf"),
  temporal_ordering_traces, width = 10, height = 6
)

temporal_ordering_traces_grid <- temporal_ordering_avg %>%
  filter(ERKi %in% c(0, 1000), directed == "directed") %>%
  mutate(log2FoldChange_rescaled_mean = if_else(
    direction_max == "pos", log2FoldChange_rescaled_mean, -1 * log2FoldChange_rescaled_mean
  )) %>%
  distinct(gene_id, Time, induction = mid_induction, ERKi, log2FoldChange_rescaled_mean) %>%
  # gather("induction", "induction_time", mid_induction, max_induction) %>%
  inner_join(
    dose_clusters %>%
      select(gene_id, dose_class = k_medoids) %>%
      filter(dose_class != "no_response_0") %>%
      mutate(dose_class = str_replace(dose_class, "_[+-]", ""))
  ) %>%
  mutate_at(vars(induction), . %>% as.character() %>% fct_reorder(as.numeric(.))) %>%
  group_nest(ERKi) %>%
  rowwise() %>%
  mutate(
    plot = {
      plot_cluster_trajectories(
        data,
        Time, log2FoldChange_rescaled_mean, Time, gene_id,
        facet_x = induction, facet_y = dose_class, all_traces = TRUE
      ) +
        labs(y = "Expression induction") +
        facet_grid(vars(dose_class), vars(induction))
    } %>%
      list()
  )

temporal_ordering_traces_grid <- temporal_ordering_avg %>%
  filter(ERKi %in% c(0, 1000), directed == "directed") %>%
  mutate(log2FoldChange = if_else(
    direction_max == "pos", log2FoldChange_rescaled_mean, -1 * log2FoldChange_rescaled_mean
  )) %>%
  distinct(gene_id, Time, induction = mid_induction, ERKi, log2FoldChange) %>%
  # gather("induction", "induction_time", mid_induction, max_induction) %>%
  inner_join(
    dose_clusters %>%
      select(gene_id, dose_class = k_medoids) %>%
      filter(dose_class != "no_response_0") %>%
      mutate(dose_class = str_replace(dose_class, "_[+-]", ""))
  ) %>%
  mutate_at(vars(induction), . %>% as.character() %>% fct_reorder(as.numeric(.))) %>%
  group_nest(ERKi) %>%
  rowwise() %>%
  mutate(
    plot = {
      plot_cluster_trajectories(
        data,
        Time, log2FoldChange, Time, gene_id,
        facet_x = induction, facet_y = dose_class, all_traces = TRUE
      ) +
        labs(y = "Expression induction") +
        facet_grid(vars(dose_class), vars(induction))
    } %>%
      list()
  )

x <- egg::ggarrange(
  plots = temporal_ordering_traces_grid[["plot"]],
  labels = temporal_ordering_traces_grid[["ERKi"]]
)

cowplot::ggsave2(
  file.path(wd, "time_dose_clusters_traces_scaled.pdf"),
  x, width = 6, height = 10
)

temporal_ordering_traces_grid <- temporal_ordering_avg %>%
  filter(ERKi %in% c(0, 1000), directed == "directed") %>%
  mutate(log2FoldChange = if_else(
    direction_max == "pos", log2FoldChange_rescaled_mean, -1 * log2FoldChange_rescaled_mean
  )) %>%
  distinct(gene_id, Time, induction = mid_induction, ERKi, log2FoldChange) %>%
  # gather("induction", "induction_time", mid_induction, max_induction) %>%
  inner_join(
    dose_clusters %>%
      select(gene_id, dose_class = k_medoids) %>%
      filter(dose_class != "no_response_0") %>%
      mutate(dose_class = str_replace(dose_class, "_[+-]", ""))
  ) %>%
  semi_join(
    padj_long %>%
      filter(padj < 0.05)
  ) %>%
  mutate_at(vars(induction), . %>% as.character() %>% fct_reorder(as.numeric(.))) %>%
  group_nest(ERKi) %>%
  rowwise() %>%
  mutate(
    plot = {
      plot_cluster_trajectories(
        data,
        Time, log2FoldChange, Time, gene_id,
        facet_x = induction, facet_y = dose_class, all_traces = TRUE
      ) +
        labs(y = "Expression induction") +
        facet_grid(vars(dose_class), vars(induction))
    } %>%
      list()
  )

x <- egg::ggarrange(
  plots = temporal_ordering_traces_grid[["plot"]],
  labels = temporal_ordering_traces_grid[["ERKi"]]
)

temporal_ordering_traces_grid <- temporal_ordering_avg %>%
  filter(ERKi %in% c(0, 1000), directed == "directed") %>%
  mutate(log2FoldChange = if_else(
    direction_max == "pos", log2FoldChange_rescaled_mean, -1 * log2FoldChange_rescaled_mean
  )) %>%
  distinct(gene_id, Time, induction = mid_induction, ERKi, log2FoldChange) %>%
  # gather("induction", "induction_time", mid_induction, max_induction) %>%
  inner_join(
    dose_clusters %>%
      select(gene_id, dose_class = k_medoids) %>%
      filter(dose_class != "no_response_0") %>%
      mutate(dose_class = str_replace(dose_class, "_[+-]", ""))
  ) %>%
  semi_join(
    padj_long %>%
      filter(padj < 0.05)
  ) %>%
  mutate_at(vars(induction), . %>% as.character() %>% fct_reorder(as.numeric(.))) %>%
  group_nest(ERKi) %>%
  rowwise() %>%
  mutate(
    plot = {
      plot_cluster_trajectories(
        data,
        Time, log2FoldChange, Time, gene_id,
        facet_x = induction, facet_y = dose_class, all_traces = TRUE
      ) +
        labs(y = "Expression induction") +
        facet_grid(vars(dose_class), vars(induction))
    } %>%
      list()
  )

# Ordered heatmaps of gene induction -------------------------------------------
###############################################################################T

# lfc_long_clusters <- lfc_long %>%
#   inner_join(condition_meta, by = "condition") %>%
#   mutate_at(vars(ERKi), as.character) %>%
#   inner_join(
#     temporal_ordering %>%
#       # rename(induction_ERKi = ERKi) %>%
#       filter(ERKi %in% c("0", "1000")),
#     by = c("gene_id", "gene_name", "ERKi")
#   )

plot_induction_heatmap <- function(
  df, row_id, row_order, col_id, col_id_order, fill, ...
) {
# plot_induction_heatmap <- function(df, row_id, col_id, fill) {
  row_id_quo <- enquo(row_id)
  row_order_quo <- enquo(row_order)
  col_id_quo <- enquo(col_id)
  col_order_quo <- enquo(col_id_order)
  fill_quo <- enquo(fill)
  row_id_order <- df %>%
    arrange(!!!eval(rlang::get_expr(row_order_quo), envir = df)) %>%
    pull(!!row_id_quo) %>%
    unique()
  col_id_order <- df %>%
    arrange(!!!eval(rlang::get_expr(col_order_quo), envir = df)) %>%
    pull(!!col_id_quo) %>%
    unique()
  df <- df %>%
    mutate(
      !!row_id_quo := factor(!!row_id_quo, levels = row_id_order),
      !!col_id_quo := factor(!!col_id_quo, levels = col_id_order)
    )
  abs_max <- df %>%
    pull(!!fill_quo) %>%
    abs() %>%
    quantile(.90, names = FALSE, na.rm = TRUE) %>%
    round(1)
  if (abs_max < 0.01) abs_max <- 0.01
  breaks <- seq(-abs_max, abs_max, by = signif(abs_max/10, 1))
  cmap <- rev(colorRampPalette(brewer.pal(7, "RdBu"))(length(breaks) - 1))
  # message(breaks)
  mat <- df %>%
    select(!!row_id_quo, !!col_id_quo, !!fill_quo) %>%
    spread(!!col_id_quo, !!fill_quo) %>%
    column_to_rownames(rlang::as_name(row_id_quo)) %>%
    as.matrix()
  # browser()
  pheatmap(
    mat,
    color = cmap,
    breaks = breaks,
    cluster_cols = FALSE,
    cluster_rows = FALSE,
    draw = FALSE,
    ...
  )
}



temporal_ordering_heatmap_data <- temporal_ordering_stats %>%
  mutate_at(vars(ERKi), as.character) %>%
  inner_join(temporal_ordering, by = c("gene_id", "gene_name", "ERKi", "directed")) %>%
  distinct(
    gene_id,
    directed,
    direction_max = direction_max.x,
    log2FoldChange,
    log2FoldChange_norm = log2FoldChange_rescaled * log2FoldChange_var_erki * if_else(direction_max.x == "pos", 1, -1),
    ERKi,
    Time,
    mid_induction,
    max_induction,
    variance
  )

temporal_ordering_heatmap_row_ann <- temporal_ordering %>%
  select(gene_id, ERKi, directed, max_induction, mid_induction) %>%
  mutate_at(
    vars(ends_with("induction")),
    ~factor(as.character(.x), levels = as.character(sort(unique(.x))))
  ) %>%
  inner_join(select(function_clusters, gene_id, class_combined), by = c("gene_id")) %>%
  group_nest(ERKi, directed, .key = "ann") %>%
  mutate(
    ann = map(
      ann,
      ~.x %>%
        column_to_rownames("gene_id") %>%
        as.data.frame()
    )
  )


temporal_ordering_heatmaps <- temporal_ordering_heatmap_data %>%
  group_nest(ERKi, directed) %>%
  inner_join(temporal_ordering_heatmap_row_ann, by = c("ERKi", "directed")) %>%
  mutate(
    data = map2(
      data, ann,
      ~plot_induction_heatmap(
        .x, gene_id, list(direction_max, mid_induction, desc(variance)), Time, list(Time), log2FoldChange,
        show_rownames = FALSE, annotation_row = .y
      )
    )
  )

pwalk(
  temporal_ordering_heatmaps,
  function(ERKi, directed, data, ...) {
    ggsave(
      file.path(wd, paste0("temporal_ordering_heatmap_erki_", ERKi, "_", directed, ".pdf")),
      data,
      width = 4, height = 5
    )
  }
)

# Scatter plot of gene ranks at different ERK concentrations -------------------
###############################################################################T

scatter_mid_max <- temporal_ordering %>%
  arrange(ERKi, mid_induction, max_induction, desc(variance)) %>%
  filter(ERKi %in% c("0", "1000")) %>%
  mutate(ERKi = paste0("rank_ERKi_", ERKi)) %>%
  group_by(ERKi, directed) %>%
  mutate(rank = 1:n()) %>%
  ungroup() %>%
  group_by(gene_name, directed) %>%
  summarize(
    rank = list(set_names(rank, ERKi)),
    variance = prod(variance)
  ) %>%
  ungroup() %>%
  inner_join(function_clusters, by = "gene_name") %>%
  unnest_wider(rank) %>%
  ggplot(aes(rank_ERKi_0, rank_ERKi_1000, alpha = variance, color = class_combined)) +
    geom_point() +
    scale_alpha_continuous(trans = "log10") +
    facet_grid(vars(directed), vars(class_combined))

ggsave(
  file.path(wd, "temporal_ordering_rank_erki_0_vs_erki_1000.pdf"),
  scatter_mid_max, width = 30, height = 8
)

