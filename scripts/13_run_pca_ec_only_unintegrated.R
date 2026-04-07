#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_rds <- "results/objects/12_ec_sct_merged.rds"
  output_rds <- "results/objects/13_ec_sct_pca_unintegrated.rds"
  features_csv <- "results/pca/13_ec_pca_features_used.csv"
  variance_csv <- "results/pca/13_ec_pca_variance_summary.csv"

  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/pca", recursive = TRUE, showWarnings = FALSE)
  dir.create("plots/pca", recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 13_run_pca_ec_only_unintegrated")
  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  DefaultAssay(seu) <- "SCT"

  message_ts("Reconstructing sample-level SCT objects from merged object")
  seu.list <- SplitObject(seu, split.by = "orig.ident")

  features <- SelectIntegrationFeatures(
    object.list = seu.list,
    nfeatures = 3000
  )

  message_ts("Selected", length(features), "shared variable features for PCA")

  VariableFeatures(seu) <- features

  seu <- RunPCA(
    object = seu,
    assay = "SCT",
    features = features,
    npcs = 30,
    verbose = TRUE
  )

  saveRDS(seu, output_rds)
  message_ts("Saved PCA object:", output_rds)

  write_csv(data.frame(feature = features), features_csv)

  stdev <- seu[["pca"]]@stdev
  variance <- stdev^2
  pct_variance <- 100 * variance / sum(variance)
  cum_pct_variance <- cumsum(pct_variance)

  variance_df <- data.frame(
    PC = seq_along(stdev),
    stdev = stdev,
    variance = variance,
    pct_variance = pct_variance,
    cum_pct_variance = cum_pct_variance
  )

  write_csv(variance_df, variance_csv)
  message_ts("Saved PCA variance summary")

  p1 <- ElbowPlot(seu, ndims = 30) + ggtitle("EC-only unintegrated PCA elbow plot")
  ggsave(
    filename = "plots/pca/13_ec_elbowplot_unintegrated.png",
    plot = p1, width = 7, height = 5, dpi = 300
  )

  p2 <- ggplot(variance_df, aes(x = PC, y = pct_variance)) +
    geom_line() + geom_point() +
    labs(title = "EC-only PCA variance explained", y = "% variance") +
    theme_classic()
  ggsave(
    filename = "plots/pca/13_ec_pca_pct_variance.png",
    plot = p2, width = 7, height = 5, dpi = 300
  )

  p3 <- ggplot(variance_df, aes(x = PC, y = cum_pct_variance)) +
    geom_line() + geom_point() +
    labs(title = "EC-only PCA cumulative variance", y = "Cumulative % variance") +
    theme_classic()
  ggsave(
    filename = "plots/pca/13_ec_pca_cumulative_variance.png",
    plot = p3, width = 7, height = 5, dpi = 300
  )

  message_ts("Stage 13_run_pca_ec_only_unintegrated complete")
}

main()
