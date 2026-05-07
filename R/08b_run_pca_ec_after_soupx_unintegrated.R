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
  input_rds <- "results/objects/08a_ec_after_soupx_sct_merged.rds"
  output_rds <- "results/objects/08b_ec_after_soupx_sct_pca_unintegrated.rds"
  features_csv <- "results/pca/08b_ec_after_soupx_pca_features_used.csv"
  variance_csv <- "results/pca/08b_ec_after_soupx_pca_variance_summary.csv"

  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/pca", recursive = TRUE, showWarnings = FALSE)
  dir.create("plots/pca", recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 08b_run_pca_ec_after_soupx_unintegrated")

  seu <- readRDS(input_rds)
  message_ts("Loaded merged object:", input_rds)

  DefaultAssay(seu) <- "SCT"

  sample_ids <- sort(unique(as.character(seu$orig.ident)))
  message_ts("Samples detected:", paste(sample_ids, collapse = ", "))

  sct_paths <- file.path("results/sct", paste0("08a_", sample_ids, "_EC_after_SoupX_SCT.rds"))
  missing_paths <- sct_paths[!file.exists(sct_paths)]
  if (length(missing_paths) > 0) {
    stop(
      paste(
        "Missing per-sample SCT objects:",
        paste(missing_paths, collapse = ", ")
      )
    )
  }

  message_ts("Loading per-sample SCT objects from results/sct")
  seu.list <- lapply(sct_paths, readRDS)
  names(seu.list) <- sample_ids

  message_ts("Selecting shared variable features from per-sample SCT objects")
  features <- SelectIntegrationFeatures(
    object.list = seu.list,
    nfeatures = 3000
  )

  message_ts("Selected", length(features), "shared variable features for PCA")

  VariableFeatures(seu) <- features

  message_ts("Running PCA on merged SCT object")
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

  p1 <- ElbowPlot(seu, ndims = 30) + ggtitle("EC after SoupX unintegrated PCA elbow plot")
  ggsave(
    filename = "plots/pca/08b_ec_after_soupx_elbowplot_unintegrated.png",
    plot = p1,
    width = 7,
    height = 5,
    dpi = 300
  )

  p2 <- ggplot(variance_df, aes(x = PC, y = pct_variance)) +
    geom_line() +
    geom_point() +
    labs(title = "EC after SoupX PCA variance explained", y = "% variance") +
    theme_classic()

  ggsave(
    filename = "plots/pca/08b_ec_after_soupx_pca_pct_variance.png",
    plot = p2,
    width = 7,
    height = 5,
    dpi = 300
  )

  p3 <- ggplot(variance_df, aes(x = PC, y = cum_pct_variance)) +
    geom_line() +
    geom_point() +
    labs(title = "EC after SoupX PCA cumulative variance", y = "Cumulative % variance") +
    theme_classic()

  ggsave(
    filename = "plots/pca/08b_ec_after_soupx_pca_cumulative_variance.png",
    plot = p3,
    width = 7,
    height = 5,
    dpi = 300
  )

  message_ts("Stage 08b_run_pca_ec_after_soupx_unintegrated complete")
}

main()
