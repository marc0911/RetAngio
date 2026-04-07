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
  input_rds <- "results/objects/17_ec_integrated_clustered_umap.rds"
  output_all <- "results/markers/18_ec_integrated_markers_all_clusters.csv"
  output_top10 <- "results/markers/18_ec_integrated_markers_top10_per_cluster.csv"
  output_sizes <- "results/markers/18_ec_integrated_cluster_sizes.csv"

  dir.create("results/markers", recursive = TRUE, showWarnings = FALSE)

  future::plan("sequential")
  options(future.globals.maxSize = 8 * 1024^3)

  message_ts("Starting stage 18_markers_ec_only_integrated_RNA")
  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  DefaultAssay(seu) <- "RNA"

  if ("RNA" %in% names(seu@assays)) {
    message_ts("Joining RNA layers before marker analysis")
    seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  }

  seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)

  Idents(seu) <- seu$SCT_snn_res.0.4

  cluster_sizes <- as.data.frame(table(Idents(seu)), stringsAsFactors = FALSE)
  colnames(cluster_sizes) <- c("cluster", "n_cells")
  cluster_sizes$cluster <- as.character(cluster_sizes$cluster)
  write_csv(cluster_sizes, output_sizes)

  markers <- FindAllMarkers(
    object = seu,
    assay = "RNA",
    only.pos = TRUE,
    test.use = "wilcox",
    logfc.threshold = 0.25,
    min.pct = 0.10,
    return.thresh = 0.05,
    verbose = TRUE
  )

  write_csv(markers, output_all)

  top10 <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
    ungroup()

  write_csv(top10, output_top10)

  message_ts("Saved all markers:", output_all)
  message_ts("Saved top10 markers:", output_top10)
  message_ts("Stage 18_markers_ec_only_integrated_RNA complete")
}

main()
