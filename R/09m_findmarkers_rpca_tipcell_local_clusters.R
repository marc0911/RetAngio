suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(readr)
  library(future)
})

log_message <- function(...) {
  msg <- paste0(...)
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg, "\n", sep = "")
}

main <- function() {
  log_message("Starting stage 09m_findmarkers_rpca_tipcell_local_clusters")

  future::plan("sequential")
  options(future.globals.maxSize = 8 * 1024^3)

  input_rds <- "results/objects/09k_rpca_tipcell_resolution_sweep.rds"
  output_dir_markers <- "results/markers"
  output_dir_review <- "results/review"

  dir.create(output_dir_markers, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir_review, recursive = TRUE, showWarnings = FALSE)

  cluster_col <- "tip_res_0_20"

  seu <- readRDS(input_rds)
  log_message("Loaded object: ", input_rds)

  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop("Cluster column not found in metadata: ", cluster_col)
  }

  clust <- seu@meta.data[[cluster_col]]
  if (all(is.na(clust))) {
    stop("All values are NA in cluster column: ", cluster_col)
  }

  clust_chr <- as.character(clust)
  clust_levels <- sort(unique(clust_chr[!is.na(clust_chr)]))
  seu@meta.data[[cluster_col]] <- factor(clust_chr, levels = clust_levels)

  Idents(seu) <- seu@meta.data[[cluster_col]]
  log_message("Using cluster column: ", cluster_col)
  log_message("Cluster levels: ", paste(levels(Idents(seu)), collapse = ", "))

  if (!"RNA" %in% Assays(seu)) {
    stop("RNA assay not found in object")
  }

  DefaultAssay(seu) <- "RNA"

  log_message("Joining RNA layers")
  seu <- JoinLayers(seu, assay = "RNA")

  log_message("Normalizing joined RNA assay")
  seu <- NormalizeData(
    object = seu,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  log_message("Running FindAllMarkers on RNA assay")
  markers <- FindAllMarkers(
    object = seu,
    assay = "RNA",
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    test.use = "wilcox",
    verbose = FALSE
  )

  markers <- markers %>%
    arrange(cluster, desc(avg_log2FC), p_val_adj, desc(pct.1))

  out_full <- file.path(output_dir_markers, "09m_rpca_tipcell_res_0_20_findallmarkers.csv")
  write_csv(markers, out_full)
  log_message("Saved full markers table: ", out_full)

  top10 <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
    ungroup()

  out_top10 <- file.path(output_dir_markers, "09m_rpca_tipcell_res_0_20_top10_markers.csv")
  write_csv(top10, out_top10)
  log_message("Saved top-10 markers table: ", out_top10)

  top20 <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>%
    ungroup()

  out_top20 <- file.path(output_dir_review, "09m_rpca_tipcell_res_0_20_top20_markers.csv")
  write_csv(top20, out_top20)
  log_message("Saved top-20 markers table: ", out_top20)

  cluster_sizes <- seu@meta.data %>%
    dplyr::count(.data[[cluster_col]], name = "n_cells") %>%
    dplyr::rename(cluster = 1)

  out_sizes <- file.path(output_dir_review, "09m_rpca_tipcell_res_0_20_cluster_sizes.csv")
  write_csv(cluster_sizes, out_sizes)
  log_message("Saved cluster sizes table: ", out_sizes)

  log_message("Number of marker rows: ", nrow(markers))
  log_message("Number of clusters with markers: ", dplyr::n_distinct(markers$cluster))
  log_message("Stage 09m_findmarkers_rpca_tipcell_local_clusters complete")
}

main()
