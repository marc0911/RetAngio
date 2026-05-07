#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

repo_root <- "/work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo"
setwd(repo_root)

input_rds <- "results/objects/05_soupx_corrected_merged.rds"
output_rds <- "results/objects/06_soupx_unintegrated_clustered_umap.rds"

dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
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

log_message("Loading SoupX-corrected merged object")
stop_if_missing(input_rds, "SoupX-corrected merged object")
obj <- readRDS(input_rds)

if (!"RNA" %in% Assays(obj)) {
  stop("RNA assay not found in input object.", call. = FALSE)
}

DefaultAssay(obj) <- "RNA"

log_message("Normalizing RNA assay")
obj <- NormalizeData(
  object = obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

log_message("Finding variable features")
obj <- FindVariableFeatures(
  object = obj,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = TRUE
)

log_message("Scaling data")
obj <- ScaleData(
  object = obj,
  features = VariableFeatures(obj),
  verbose = TRUE
)

log_message("Running PCA")
obj <- RunPCA(
  object = obj,
  features = VariableFeatures(obj),
  npcs = 30,
  verbose = TRUE
)

log_message("Finding neighbors")
obj <- FindNeighbors(
  object = obj,
  dims = 1:20,
  verbose = TRUE
)

log_message("Finding clusters")
obj <- FindClusters(
  object = obj,
  resolution = 0.4,
  verbose = TRUE
)

log_message("Running UMAP")
obj <- RunUMAP(
  object = obj,
  dims = 1:20,
  reduction = "pca",
  reduction.name = "umap.unintegrated",
  reduction.key = "UMAPUNINT_",
  verbose = TRUE
)

obj$stage06_cluster <- as.character(Idents(obj))

if (is.null(obj@misc)) {
  obj@misc <- list()
}

obj@misc$stage06_unintegrated_clustering <- list(
  source_object = input_rds,
  assay = "RNA",
  normalization = "LogNormalize",
  variable_features_n = 3000,
  pca_npcs = 30,
  neighbor_dims = 20,
  clustering_resolution = 0.4,
  umap_reduction_name = "umap.unintegrated",
  run_datetime = as.character(Sys.time())
)

saveRDS(obj, output_rds)

log_message("Saved unintegrated clustered object: ", output_rds)
log_message("Number of cells: ", ncol(obj))
log_message("Number of clusters: ", length(unique(Idents(obj))))
log_message("Cluster levels: ", paste(levels(Idents(obj)), collapse = ", "))
log_message("Stage 06 completed successfully")
