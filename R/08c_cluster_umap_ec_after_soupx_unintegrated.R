#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_rds <- "results/objects/08b_ec_after_soupx_sct_pca_unintegrated.rds"
  decision_file <- "decisions/08b_ec_after_soupx_unintegrated_pca_dims.yml"
  output_rds <- "results/objects/08c_ec_after_soupx_unintegrated_clustered_umap.rds"

  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 08c_cluster_umap_ec_after_soupx_unintegrated")

  dec <- yaml::read_yaml(decision_file)
  if (is.null(dec$selected_dims)) {
    stop("Decision file does not contain 'selected_dims'.")
  }

  dims_use <- as.integer(unlist(dec$selected_dims))
  if (length(dims_use) == 0) {
    stop("No dimensions found in decision file.")
  }

  message_ts("Using dims from decision file:", paste(dims_use, collapse = ","))

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  DefaultAssay(seu) <- "SCT"

  message_ts("Running FindNeighbors")
  seu <- FindNeighbors(
    object = seu,
    reduction = "pca",
    dims = dims_use,
    verbose = TRUE
  )

  message_ts("Running FindClusters at resolution 0.4")
  seu <- FindClusters(
    object = seu,
    resolution = 0.4,
    verbose = TRUE
  )

  message_ts("Running UMAP")
  seu <- RunUMAP(
    object = seu,
    reduction = "pca",
    dims = dims_use,
    reduction.name = "umap.unintegrated",
    reduction.key = "UMAPUNINT_",
    verbose = TRUE
  )

  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    seu$seurat_clusters <- as.character(Idents(seu))
  }

  seu$stage08c_cluster <- as.character(Idents(seu))

  if (is.null(seu@misc)) {
    seu@misc <- list()
  }

  seu@misc$stage08c_cluster_umap_unintegrated <- list(
    source_object = input_rds,
    decision_file = decision_file,
    dims_used = dims_use,
    clustering_resolution = 0.4,
    umap_reduction_name = "umap.unintegrated",
    run_datetime = as.character(Sys.time())
  )

  saveRDS(seu, output_rds)

  message_ts("Saved clustered UMAP object:", output_rds)
  message_ts("Number of cells:", ncol(seu))
  message_ts("Number of clusters:", length(unique(Idents(seu))))
  message_ts("Cluster levels:", paste(levels(Idents(seu)), collapse = ", "))
  message_ts("Stage 08c_cluster_umap_ec_after_soupx_unintegrated complete")
}

main()
