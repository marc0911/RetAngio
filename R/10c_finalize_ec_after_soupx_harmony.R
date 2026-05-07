#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  decision_file <- "decisions/10c_ec_after_soupx_harmony_resolution_choice.yaml"
  output_rds <- "results/objects/10c_ec_after_soupx_harmony_final.rds"

  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 10c_finalize_ec_after_soupx_harmony")
  message_ts("Reading decision file")

  dec <- yaml::read_yaml(decision_file)

  input_rds <- dec$source_object
  chosen_resolution <- as.character(dec$chosen_resolution)
  chosen_cluster_column <- dec$chosen_cluster_column
  chosen_reduction <- dec$chosen_reduction
  chosen_umap <- dec$chosen_umap

  message_ts("Using resolution:", chosen_resolution)
  message_ts("Using cluster column:", chosen_cluster_column)

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  if (!chosen_cluster_column %in% colnames(seu@meta.data)) {
    stop(paste("Chosen cluster column not found:", chosen_cluster_column))
  }

  if (!chosen_reduction %in% names(seu@reductions)) {
    stop(paste("Chosen reduction not found:", chosen_reduction))
  }

  if (!chosen_umap %in% names(seu@reductions)) {
    stop(paste("Chosen UMAP reduction not found:", chosen_umap))
  }

  seu$seurat_clusters <- as.character(seu@meta.data[[chosen_cluster_column]])
  cluster_levels <- as.character(sort(unique(as.numeric(seu$seurat_clusters))))
  seu$seurat_clusters <- factor(seu$seurat_clusters, levels = cluster_levels)
  Idents(seu) <- "seurat_clusters"

  if (is.null(seu@misc)) {
    seu@misc <- list()
  }

  seu@misc$stage10c_harmony_finalize <- list(
    decision_file = decision_file,
    source_object = input_rds,
    chosen_resolution = chosen_resolution,
    chosen_cluster_column = chosen_cluster_column,
    chosen_reduction = chosen_reduction,
    chosen_umap = chosen_umap,
    cluster_levels = cluster_levels,
    run_datetime = as.character(Sys.time())
  )

  saveRDS(seu, output_rds)

  message_ts("Saved final Harmony object:", output_rds)
  message_ts("Number of cells:", ncol(seu))
  message_ts("Number of clusters:", length(cluster_levels))
  message_ts("Cluster levels:", paste(cluster_levels, collapse = ", "))
  message_ts("Stage 10c_finalize_ec_after_soupx_harmony complete")
}

main()
