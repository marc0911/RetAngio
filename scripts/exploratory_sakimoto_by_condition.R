#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_rds <- "results/objects/14_ec_unintegrated_clustered_umap.rds"

  plot_dir <- "plots/exploratory_sakimoto_by_condition"
  result_dir <- "results/exploratory_sakimoto_by_condition"

  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting exploratory Sakimoto distribution by condition")
  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  DefaultAssay(seu) <- "RNA"
  seu <- JoinLayers(seu, assay = "RNA")
  seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)

  genes_needed <- c("Bst1", "Procr")
  missing_genes <- setdiff(genes_needed, rownames(seu))
  if (length(missing_genes) > 0) {
    stop("Missing required genes: ", paste(missing_genes, collapse = ", "))
  }

  cluster_col <- "SCT_snn_res.0.4"
  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop("Expected cluster column not found: ", cluster_col)
  }

  df <- FetchData(seu, vars = c("orig.ident", cluster_col, "Bst1", "Procr")) %>%
    rename(cluster = all_of(cluster_col)) %>%
    mutate(
      Bst1_positive = Bst1 > 0,
      Procr_positive = Procr > 0,
      double_positive = Bst1_positive & Procr_positive,
      Procr_only = Procr_positive & !Bst1_positive,
      Bst1_only = Bst1_positive & !Procr_positive
    )

  write_csv(
    df,
    file.path(result_dir, "sakimoto_condition_cell_level_table.csv")
  )

  summary_by_condition <- df %>%
    group_by(orig.ident) %>%
    summarise(
      n_cells = n(),
      n_Bst1_positive = sum(Bst1_positive),
      pct_Bst1_positive = 100 * mean(Bst1_positive),
      n_Procr_positive = sum(Procr_positive),
      pct_Procr_positive = 100 * mean(Procr_positive),
      n_double_positive = sum(double_positive),
      pct_double_positive = 100 * mean(double_positive),
      n_Procr_only = sum(Procr_only),
      pct_Procr_only = 100 * mean(Procr_only),
      n_Bst1_only = sum(Bst1_only),
      pct_Bst1_only = 100 * mean(Bst1_only),
      avg_Bst1 = mean(Bst1),
      avg_Procr = mean(Procr),
      .groups = "drop"
    )

  write_csv(
    summary_by_condition,
    file.path(result_dir, "sakimoto_distribution_by_condition.csv")
  )

  summary_by_condition_cluster <- df %>%
    group_by(orig.ident, cluster) %>%
    summarise(
      n_cells = n(),
      n_Bst1_positive = sum(Bst1_positive),
      pct_Bst1_positive = 100 * mean(Bst1_positive),
      n_Procr_positive = sum(Procr_positive),
      pct_Procr_positive = 100 * mean(Procr_positive),
      n_double_positive = sum(double_positive),
      pct_double_positive = 100 * mean(double_positive),
      n_Procr_only = sum(Procr_only),
      pct_Procr_only = 100 * mean(Procr_only),
      n_Bst1_only = sum(Bst1_only),
      pct_Bst1_only = 100 * mean(Bst1_only),
      avg_Bst1 = mean(Bst1),
      avg_Procr = mean(Procr),
      .groups = "drop"
    ) %>%
    arrange(orig.ident, suppressWarnings(as.numeric(as.character(cluster))), cluster)

  write_csv(
    summary_by_condition_cluster,
    file.path(result_dir, "sakimoto_distribution_by_condition_and_cluster.csv")
  )

  p1 <- ggplot(summary_by_condition, aes(x = orig.ident, y = pct_Procr_positive)) +
    geom_col() +
    theme_bw() +
    labs(
      title = "Procr+ ECs by condition",
      x = "Condition",
      y = "% Procr+ cells"
    )

  ggsave(
    filename = file.path(plot_dir, "sakimoto_pct_Procr_positive_by_condition.png"),
    plot = p1, width = 8, height = 6, dpi = 300
  )

  p2 <- ggplot(summary_by_condition, aes(x = orig.ident, y = pct_Bst1_positive)) +
    geom_col() +
    theme_bw() +
    labs(
      title = "Bst1+ ECs by condition",
      x = "Condition",
      y = "% Bst1+ cells"
    )

  ggsave(
    filename = file.path(plot_dir, "sakimoto_pct_Bst1_positive_by_condition.png"),
    plot = p2, width = 8, height = 6, dpi = 300
  )

  p3 <- ggplot(summary_by_condition, aes(x = orig.ident, y = pct_double_positive)) +
    geom_col() +
    theme_bw() +
    labs(
      title = "Bst1+ / Procr+ ECs by condition",
      x = "Condition",
      y = "% double-positive cells"
    )

  ggsave(
    filename = file.path(plot_dir, "sakimoto_pct_double_positive_by_condition.png"),
    plot = p3, width = 8, height = 6, dpi = 300
  )

  p4 <- ggplot(summary_by_condition, aes(x = orig.ident, y = pct_Procr_only)) +
    geom_col() +
    theme_bw() +
    labs(
      title = "Procr+ / Bst1- ECs by condition",
      x = "Condition",
      y = "% Procr+ / Bst1- cells"
    )

  ggsave(
    filename = file.path(plot_dir, "sakimoto_pct_Procr_only_by_condition.png"),
    plot = p4, width = 8, height = 6, dpi = 300
  )

  p5 <- ggplot(summary_by_condition_cluster, aes(x = cluster, y = pct_double_positive, fill = orig.ident)) +
    geom_col(position = "dodge") +
    theme_bw() +
    labs(
      title = "Bst1+ / Procr+ ECs by condition and cluster",
      x = "Cluster",
      y = "% double-positive cells"
    )

  ggsave(
    filename = file.path(plot_dir, "sakimoto_pct_double_positive_by_condition_and_cluster.png"),
    plot = p5, width = 12, height = 6, dpi = 300
  )

  p6 <- ggplot(summary_by_condition_cluster, aes(x = cluster, y = pct_Procr_only, fill = orig.ident)) +
    geom_col(position = "dodge") +
    theme_bw() +
    labs(
      title = "Procr+ / Bst1- ECs by condition and cluster",
      x = "Cluster",
      y = "% Procr+ / Bst1- cells"
    )

  ggsave(
    filename = file.path(plot_dir, "sakimoto_pct_Procr_only_by_condition_and_cluster.png"),
    plot = p6, width = 12, height = 6, dpi = 300
  )

  message_ts("Exploratory Sakimoto distribution by condition complete")
}

main()
