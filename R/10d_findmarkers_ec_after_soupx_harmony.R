#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_rds <- "results/objects/10c_ec_after_soupx_harmony_final.rds"
  output_markers <- "results/markers/10d_ec_after_soupx_harmony_findallmarkers.csv"
  output_top10 <- "results/markers/10d_ec_after_soupx_harmony_top10_markers.csv"

  dir.create("results/markers", recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 10d_findmarkers_ec_after_soupx_harmony")

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  if (!"RNA" %in% Assays(seu)) {
    stop("RNA assay not found in object.")
  }

  DefaultAssay(seu) <- "RNA"

  message_ts("Joining RNA layers")
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])

  message_ts("Normalizing joined RNA assay")
  seu <- NormalizeData(
    object = seu,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = TRUE
  )

  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    seu$seurat_clusters <- as.character(Idents(seu))
  }
  Idents(seu) <- "seurat_clusters"

  message_ts("Running FindAllMarkers on RNA assay")
  markers <- FindAllMarkers(
    object = seu,
    assay = "RNA",
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    verbose = TRUE
  )

  write_csv(markers, output_markers)

  top10 <- markers %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup()

  write_csv(top10, output_top10)

  message_ts("Saved full markers table:", output_markers)
  message_ts("Saved top-10 markers table:", output_top10)
  message_ts("Number of marker rows:", nrow(markers))
  message_ts("Number of clusters with markers:", dplyr::n_distinct(markers$cluster))
  message_ts("Stage 10d_findmarkers_ec_after_soupx_harmony complete")
}

main()
