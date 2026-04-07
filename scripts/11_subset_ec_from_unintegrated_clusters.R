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
  input_raw <- "results/objects/04_qc_doublet_filtered.rds"
  input_clustered <- "results/objects/08_unintegrated_clustered_umap.rds"
  output_rds <- "results/objects/11_ec_enriched_raw.rds"
  output_cells <- "results/metadata/11_ec_enriched_cells.csv"
  output_summary <- "results/metadata/11_ec_enriched_summary.csv"

  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/metadata", recursive = TRUE, showWarnings = FALSE)

  keep_clusters <- c("0","1","2","3","4","5","6")

  message_ts("Starting stage 11_subset_ec_from_unintegrated_clusters")

  message_ts("Loading clustered object:", input_clustered)
  seu_clustered <- readRDS(input_clustered)

  if (!"seurat_clusters" %in% colnames(seu_clustered@meta.data)) {
    stop("seurat_clusters not found in clustered object metadata.")
  }

  clustered_meta <- seu_clustered@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(cell, orig.ident, seurat_clusters)

  clustered_meta$seurat_clusters <- as.character(clustered_meta$seurat_clusters)

  cells_keep <- clustered_meta %>%
    dplyr::filter(seurat_clusters %in% keep_clusters) %>%
    dplyr::pull(cell)

  if (length(cells_keep) == 0) {
    stop("No cells selected for EC-enriched object.")
  }

  message_ts("Cells retained:", length(cells_keep))

  write_csv(
    clustered_meta %>%
      dplyr::mutate(retain_ec = seurat_clusters %in% keep_clusters),
    output_cells
  )

  summary_df <- clustered_meta %>%
    dplyr::count(seurat_clusters, name = "n_cells") %>%
    dplyr::mutate(retain_ec = seurat_clusters %in% keep_clusters)

  write_csv(summary_df, output_summary)

  message_ts("Loading raw filtered object:", input_raw)
  seu_raw <- readRDS(input_raw)

  common_cells <- intersect(colnames(seu_raw), cells_keep)

  if (length(common_cells) != length(cells_keep)) {
    warning("Some selected cells were not found in raw object. Keeping common cells only.")
  }

  if (length(common_cells) == 0) {
    stop("No overlapping cells between clustered object and raw filtered object.")
  }

  message_ts("Overlapping cells retained:", length(common_cells))

  seu_ec <- subset(seu_raw, cells = common_cells)

  seu_ec$stage11_source_cluster <- clustered_meta$seurat_clusters[
    match(colnames(seu_ec), clustered_meta$cell)
  ]

  saveRDS(seu_ec, output_rds)
  message_ts("Saved EC-enriched object:", output_rds)

  message_ts("Stage 11_subset_ec_from_unintegrated_clusters complete")
}

main()
