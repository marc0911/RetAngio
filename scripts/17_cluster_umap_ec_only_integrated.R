#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(readr)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

save_png <- function(plot_obj, filename, width = 10, height = 8, dpi = 300) {
  ggsave(filename = filename, plot = plot_obj, width = width, height = height, dpi = dpi)
}

main <- function() {
  input_rds <- "results/objects/16_ec_rpca_integrated.rds"
  output_rds <- "results/objects/17_ec_integrated_clustered_umap.rds"

  out_total <- "results/clusters/17_ec_integrated_cluster_counts_total.csv"
  out_by_sample <- "results/clusters/17_ec_integrated_cluster_counts_by_sample.csv"

  plot_clusters <- "plots/integration/17_ec_integrated_umap_clusters_res0.4.png"
  plot_by_sample <- "plots/integration/17_ec_integrated_umap_by_sample.png"
  plot_split <- "plots/integration/17_ec_integrated_umap_clusters_split_samples.png"

  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/clusters", recursive = TRUE, showWarnings = FALSE)
  dir.create("plots/integration", recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 17_cluster_umap_ec_only_integrated")
  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  if (!"integrated.rpca" %in% names(seu@reductions)) {
    stop("integrated.rpca reduction not found")
  }

  dims_use <- 1:12
  res_grid <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)

  message_ts("Using dims:", paste(dims_use, collapse = ", "))
  message_ts("Using resolution grid:", paste(res_grid, collapse = ", "))

  seu <- FindNeighbors(
    seu,
    reduction = "integrated.rpca",
    dims = dims_use,
    verbose = TRUE
  )

  for (res in res_grid) {
    message_ts("Running FindClusters at resolution", res)
    seu <- FindClusters(seu, resolution = res, verbose = FALSE)
  }

  seu <- RunUMAP(
    seu,
    reduction = "integrated.rpca",
    dims = dims_use,
    reduction.name = "umap.integrated",
    reduction.key = "UMAPINT_",
    verbose = TRUE
  )

  Idents(seu) <- seu$SCT_snn_res.0.4

  p1 <- DimPlot(seu, reduction = "umap.integrated", group.by = "SCT_snn_res.0.4", label = TRUE) +
    ggtitle("EC-only integrated RPCA") +
    theme_bw()

  p2 <- DimPlot(seu, reduction = "umap.integrated", group.by = "orig.ident") +
    ggtitle("EC-only integrated RPCA by sample") +
    theme_bw()

  p3 <- DimPlot(seu, reduction = "umap.integrated", group.by = "SCT_snn_res.0.4", split.by = "orig.ident", label = TRUE, ncol = 3) +
    theme_bw()

  save_png(p1, plot_clusters, width = 10, height = 8)
  save_png(p2, plot_by_sample, width = 10, height = 8)
  save_png(p3, plot_split, width = 16, height = 10)

  cluster_total <- as.data.frame(table(seu$SCT_snn_res.0.4), stringsAsFactors = FALSE)
  colnames(cluster_total) <- c("cluster", "n")
  cluster_total$cluster <- as.character(cluster_total$cluster)
  write_csv(cluster_total, out_total)

  cluster_by_sample <- seu@meta.data %>%
    count(orig.ident, SCT_snn_res.0.4, name = "n") %>%
    rename(cluster = SCT_snn_res.0.4) %>%
    arrange(orig.ident, cluster)
  write_csv(cluster_by_sample, out_by_sample)

  saveRDS(seu, output_rds)
  message_ts("Saved object:", output_rds)
  message_ts("Stage 17_cluster_umap_ec_only_integrated complete")
}

main()
