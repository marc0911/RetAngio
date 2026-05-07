#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

repo_root <- "/work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo"
setwd(repo_root)

input_rds <- "results/objects/06_soupx_unintegrated_clustered_umap.rds"
output_markers_csv <- "results/markers/07a_unintegrated_after_soupx_findallmarkers.csv"
output_top_markers_csv <- "results/markers/07a_unintegrated_after_soupx_top10_markers.csv"

dir.create("results/markers", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

log_message <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", ts, paste(..., collapse = "")))
}

stop_if_missing <- function(path, label = NULL) {
  if (!file.exists(path) && !dir.exists(path)) {
    if (is.null(label)) {
      stop(sprintf("Required path does not exist: %s", path), call. = FALSE)
    } else {
      stop(sprintf("Required %s does not exist: %s", label, path), call. = FALSE)
    }
  }
}

log_message("Loading unintegrated clustered SoupX-corrected object")
stop_if_missing(input_rds, "unintegrated clustered SoupX-corrected object")
obj <- readRDS(input_rds)

if (!"RNA" %in% Assays(obj)) {
  stop("RNA assay not found in input object.", call. = FALSE)
}

DefaultAssay(obj) <- "RNA"

log_message("Joining RNA layers")
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

log_message("Normalizing joined RNA assay")
obj <- NormalizeData(
  object = obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
  obj$seurat_clusters <- as.character(Idents(obj))
}

Idents(obj) <- "seurat_clusters"

log_message("Running FindAllMarkers on RNA assay")
markers <- FindAllMarkers(
  object = obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  verbose = TRUE
)

if (nrow(markers) == 0) {
  stop("FindAllMarkers returned zero rows.", call. = FALSE)
}

markers <- markers %>%
  arrange(cluster, desc(avg_log2FC), p_val_adj)

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(cluster, desc(avg_log2FC), p_val_adj)

write.csv(markers, output_markers_csv, row.names = FALSE)
write.csv(top_markers, output_top_markers_csv, row.names = FALSE)

log_message("Saved full markers table: ", output_markers_csv)
log_message("Saved top-10 markers table: ", output_top_markers_csv)
log_message("Number of marker rows: ", nrow(markers))
log_message("Number of clusters with markers: ", length(unique(markers$cluster)))
log_message("Stage 07a completed successfully")
