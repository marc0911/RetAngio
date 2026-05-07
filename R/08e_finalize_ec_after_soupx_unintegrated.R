#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_rds <- "results/objects/08d_ec_after_soupx_unintegrated_resolution_sweep.rds"
  decision_file <- "decisions/08d_ec_after_soupx_unintegrated_clustering_resolution.yml"
  output_rds <- "results/objects/08e_ec_after_soupx_unintegrated_final.rds"

  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 08e_finalize_ec_after_soupx_unintegrated")

  dec <- yaml::read_yaml(decision_file)

  if (is.null(dec$selected_parameters$clustering_resolution)) {
    stop("Decision file missing selected_parameters$clustering_resolution")
  }
  if (is.null(dec$selected_parameters$dims_used)) {
    stop("Decision file missing selected_parameters$dims_used")
  }

  res_use <- as.numeric(dec$selected_parameters$clustering_resolution)
  dims_use <- as.integer(unlist(dec$selected_parameters$dims_used))
  cluster_col <- paste0("SCT_snn_res.", res_use)

  message_ts("Using resolution:", res_use)
  message_ts("Using dims:", paste(dims_use, collapse = ","))

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop(paste("Cluster column not found in object:", cluster_col))
  }

  Idents(seu) <- cluster_col
  seu$seurat_clusters <- as.character(Idents(seu))
  seu$final_unintegrated_cluster <- as.character(Idents(seu))

  if (is.null(seu@misc)) {
    seu@misc <- list()
  }

  seu@misc$stage08e_finalize_unintegrated <- list(
    source_object = input_rds,
    decision_file = decision_file,
    dims_used = dims_use,
    clustering_resolution = res_use,
    cluster_column = cluster_col,
    reduction_used = "umap.unintegrated",
    run_datetime = as.character(Sys.time())
  )

  saveRDS(seu, output_rds)

  message_ts("Saved final unintegrated object:", output_rds)
  message_ts("Number of cells:", ncol(seu))
  message_ts("Number of clusters:", length(unique(Idents(seu))))
  message_ts("Cluster levels:", paste(levels(Idents(seu)), collapse = ", "))
  message_ts("Stage 08e_finalize_ec_after_soupx_unintegrated complete")
}

main()
