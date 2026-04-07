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
  input_rds <- "results/objects/17_ec_integrated_clustered_umap.rds"

  dir.create("plots/integration", recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 17b_replot_ec_only_integrated")
  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  main_res_col <- "SCT_snn_res.0.4"
  if (!main_res_col %in% colnames(seu@meta.data)) {
    stop("Expected main resolution column not found: ", main_res_col)
  }

  if (!"umap.integrated" %in% names(seu@reductions)) {
    stop("Expected reduction not found: umap.integrated")
  }

  Idents(seu) <- seu[[main_res_col, drop = TRUE]]

  p1 <- DimPlot(
    seu,
    reduction = "umap.integrated",
    group.by = main_res_col,
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) + ggtitle("EC-only integrated UMAP (res 0.4)")

  ggsave(
    filename = "plots/integration/17_ec_integrated_umap_clusters_res0.4.png",
    plot = p1, width = 8, height = 6, dpi = 300
  )

  p2 <- DimPlot(
    seu,
    reduction = "umap.integrated",
    group.by = "orig.ident",
    raster = FALSE
  ) + ggtitle("EC-only integrated UMAP by sample")

  ggsave(
    filename = "plots/integration/17_ec_integrated_umap_by_sample.png",
    plot = p2, width = 8, height = 6, dpi = 300
  )

  p3 <- DimPlot(
    seu,
    reduction = "umap.integrated",
    group.by = main_res_col,
    split.by = "orig.ident",
    label = FALSE,
    raster = FALSE
  ) + ggtitle("EC-only integrated UMAP split by sample")

  ggsave(
    filename = "plots/integration/17_ec_integrated_umap_clusters_split_samples.png",
    plot = p3, width = 16, height = 8, dpi = 300
  )

  message_ts("Stage 17b_replot_ec_only_integrated complete")
}

main()
