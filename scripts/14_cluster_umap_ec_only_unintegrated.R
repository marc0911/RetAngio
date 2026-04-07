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
  input_rds <- "results/objects/13_ec_sct_pca_unintegrated.rds"
  output_rds <- "results/objects/14_ec_unintegrated_clustered_umap.rds"
  counts_total_csv <- "results/clusters/14_ec_unintegrated_cluster_counts_total.csv"
  counts_by_sample_csv <- "results/clusters/14_ec_unintegrated_cluster_counts_by_sample.csv"

  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/clusters", recursive = TRUE, showWarnings = FALSE)
  dir.create("plots/umap", recursive = TRUE, showWarnings = FALSE)

  dims_use <- 1:12
  res_grid <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)

  message_ts("Starting stage 14_cluster_umap_ec_only_unintegrated")
  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  DefaultAssay(seu) <- "SCT"

  message_ts("Using dims:", paste(dims_use, collapse = ", "))
  message_ts("Using resolution grid:", paste(res_grid, collapse = ", "))

  seu <- FindNeighbors(
    seu,
    reduction = "pca",
    dims = dims_use,
    verbose = TRUE
  )

  for (res in res_grid) {
    message_ts("Running FindClusters at resolution", res)
    seu <- FindClusters(
      seu,
      resolution = res,
      verbose = FALSE
    )
  }

  main_res_col <- "SCT_snn_res.0.4"
  if (!main_res_col %in% colnames(seu@meta.data)) {
    stop("Expected main resolution column not found: ", main_res_col)
  }

  Idents(seu) <- seu[[main_res_col, drop = TRUE]]

  seu <- RunUMAP(
    seu,
    reduction = "pca",
    dims = dims_use,
    reduction.name = "umap",
    reduction.key = "UMAP_",
    verbose = TRUE
  )

  saveRDS(seu, output_rds)
  message_ts("Saved object:", output_rds)

  total_df <- as.data.frame(table(seu[[main_res_col, drop = TRUE]]), stringsAsFactors = FALSE)
  colnames(total_df) <- c("cluster", "n")
  total_df$cluster <- as.character(total_df$cluster)
  write_csv(total_df, counts_total_csv)

  by_sample_df <- seu@meta.data %>%
    mutate(cluster = .data[[main_res_col]]) %>%
    count(orig.ident, cluster, name = "n") %>%
    arrange(orig.ident, cluster)
  write_csv(by_sample_df, counts_by_sample_csv)

  p1 <- DimPlot(
    seu,
    reduction = "umap",
    group.by = main_res_col,
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) + ggtitle("EC-only unintegrated UMAP (res 0.4)")

  ggsave(
    filename = "plots/umap/14_ec_umap_unintegrated_res0.4.png",
    plot = p1, width = 8, height = 6, dpi = 300
  )

  p2 <- DimPlot(
    seu,
    reduction = "umap",
    group.by = "orig.ident",
    raster = FALSE
  ) + ggtitle("EC-only unintegrated UMAP by sample")

  ggsave(
    filename = "plots/umap/14_ec_umap_unintegrated_by_sample.png",
    plot = p2, width = 8, height = 6, dpi = 300
  )

  p3 <- DimPlot(
    seu,
    reduction = "umap",
    group.by = main_res_col,
    split.by = "orig.ident",
    label = FALSE,
    raster = FALSE
  ) + ggtitle("EC-only unintegrated UMAP split by sample")

  ggsave(
    filename = "plots/umap/14_ec_umap_unintegrated_res0.4_split_samples.png",
    plot = p3, width = 16, height = 8, dpi = 300
  )

  message_ts("Stage 14_cluster_umap_ec_only_unintegrated complete")
}

main()
