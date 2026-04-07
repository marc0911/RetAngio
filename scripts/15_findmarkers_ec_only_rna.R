#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(future)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_rds <- "results/objects/14_ec_unintegrated_clustered_umap.rds"
  output_rds <- "results/objects/15_ec_unintegrated_rna_ready.rds"
  markers_all_csv <- "results/markers/15_ec_markers_all_clusters.csv"
  markers_top10_csv <- "results/markers/15_ec_markers_top10_per_cluster.csv"
  cluster_sizes_csv <- "results/markers/15_ec_cluster_sizes.csv"
  cluster_by_sample_csv <- "results/markers/15_ec_cluster_counts_by_sample.csv"

  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/markers", recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 15_findmarkers_ec_only_rna")

  future::plan("sequential")
  options(future.globals.maxSize = 8 * 1024^3)

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  main_res_col <- "SCT_snn_res.0.4"
  if (!main_res_col %in% colnames(seu@meta.data)) {
    stop("Expected clustering column not found: ", main_res_col)
  }

  Idents(seu) <- seu[[main_res_col, drop = TRUE]]

  DefaultAssay(seu) <- "RNA"
  message_ts("Default assay set to RNA")

  message_ts("Joining RNA layers")
  seu <- JoinLayers(seu)

  message_ts("Normalizing RNA assay")
  seu <- NormalizeData(seu, verbose = FALSE)

  saveRDS(seu, output_rds)
  message_ts("Saved RNA-ready object:", output_rds)

  cluster_sizes <- as.data.frame(table(Idents(seu)), stringsAsFactors = FALSE)
  colnames(cluster_sizes) <- c("cluster", "n")
  cluster_sizes$cluster <- as.character(cluster_sizes$cluster)
  write_csv(cluster_sizes, cluster_sizes_csv)

  by_sample <- seu@meta.data %>%
    mutate(cluster = .data[[main_res_col]]) %>%
    count(orig.ident, cluster, name = "n") %>%
    arrange(orig.ident, cluster)
  write_csv(by_sample, cluster_by_sample_csv)

  message_ts("Running FindAllMarkers on RNA assay")
  markers <- FindAllMarkers(
    object = seu,
    assay = "RNA",
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    verbose = TRUE
  )

  write_csv(markers, markers_all_csv)

  fc_col <- if ("avg_log2FC" %in% colnames(markers)) "avg_log2FC" else
    if ("avg_logFC" %in% colnames(markers)) "avg_logFC" else
      stop("No logFC column found in markers table")

  markers_top10 <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = .data[[fc_col]], n = 10, with_ties = FALSE) %>%
    ungroup()

  write_csv(markers_top10, markers_top10_csv)

  message_ts("Stage 15_findmarkers_ec_only_rna complete")
}

main()
