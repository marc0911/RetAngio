source("R/00_utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tibble)
  library(yaml)
})

main <- function() {
  log_message("Starting stage 08_unintegrated_clustering_umap")

  cfg <- read_config()
  outdir <- cfg$project$output_dir
  sample_col <- cfg$analysis$sample_column

  ensure_dir(file.path(outdir, "objects"))
  ensure_dir(file.path(outdir, "clusters"))
  ensure_dir(file.path("plots", "umap"))

  in_rds <- file.path(outdir, "objects", "07_sct_pca_unintegrated.rds")
  dec_path <- "decisions/05_unintegrated_dims_resolution.yml"

  if (!file.exists(in_rds)) stop("Missing input object: ", in_rds)
  if (!file.exists(dec_path)) stop("Missing decision file: ", dec_path)

  seu <- readRDS(in_rds)
  dec <- yaml::read_yaml(dec_path)

  dims_vec <- seq(
    dec$decision$unintegrated_dims_primary[[1]],
    dec$decision$unintegrated_dims_primary[[2]]
  )
  res_grid <- unlist(dec$decision$resolution_grid)

  log_message("Using dims: ", paste(range(dims_vec), collapse = ":"))
  log_message("Using resolution grid: ", paste(res_grid, collapse = ", "))

  DefaultAssay(seu) <- "SCT"

  seu <- FindNeighbors(
    seu,
    reduction = "pca",
    dims = dims_vec,
    graph.name = "unintegrated_snn",
    verbose = TRUE
  )

  for (res in res_grid) {
    log_message("Running FindClusters at resolution ", res)
    seu <- FindClusters(
      seu,
      graph.name = "unintegrated_snn",
      resolution = res,
      algorithm = 1,
      verbose = FALSE
    )
  }

  seu <- RunUMAP(
    seu,
    reduction = "pca",
    dims = dims_vec,
    reduction.name = "umap.unintegrated",
    reduction.key = "UMAPUN_",
    verbose = TRUE
  )

  # Main resolution for quick presentation
  main_res <- "unintegrated_snn_res.0.4"
  if (!main_res %in% colnames(seu@meta.data)) {
    stop(
      "Expected main resolution column not found: ", main_res,
      "\nAvailable resolution columns: ",
      paste(grep("_res\\.", colnames(seu@meta.data), value = TRUE), collapse = ", ")
    )
  }

  seu$seurat_clusters <- seu@meta.data[[main_res]]

  saveRDS(seu, file.path(outdir, "objects", "08_unintegrated_clustered_umap.rds"))
  log_message("Saved object: results/objects/08_unintegrated_clustered_umap.rds")

  cluster_tbl <- seu@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    transmute(
      cell_id = cell_id,
      sample_id = .data[[sample_col]],
      cluster = as.character(.data[[main_res]])
    ) %>%
    count(cluster, sample_id, name = "n") %>%
    arrange(as.numeric(cluster), sample_id)

  write_csv(cluster_tbl, file.path(outdir, "clusters", "08_unintegrated_cluster_counts_by_sample.csv"))

  cluster_totals <- seu@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    transmute(cluster = as.character(.data[[main_res]])) %>%
    count(cluster, name = "n") %>%
    arrange(as.numeric(cluster))

  write_csv(cluster_totals, file.path(outdir, "clusters", "08_unintegrated_cluster_counts_total.csv"))

  p1 <- DimPlot(
    seu,
    reduction = "umap.unintegrated",
    group.by = main_res,
    label = TRUE,
    repel = TRUE
  ) + ggtitle("Unintegrated UMAP (res 0.4)")

  ggsave(
    filename = file.path("plots", "umap", "08_umap_unintegrated_res0.4.png"),
    plot = p1,
    width = 8,
    height = 6,
    dpi = 300
  )

  p2 <- DimPlot(
    seu,
    reduction = "umap.unintegrated",
    group.by = sample_col
  ) + ggtitle("Unintegrated UMAP by sample")

  ggsave(
    filename = file.path("plots", "umap", "08_umap_unintegrated_by_sample.png"),
    plot = p2,
    width = 8,
    height = 6,
    dpi = 300
  )

  p3 <- DimPlot(
    seu,
    reduction = "umap.unintegrated",
    group.by = main_res,
    split.by = sample_col,
    ncol = 3,
    label = FALSE
  ) + ggtitle("Unintegrated UMAP split by sample (res 0.4)")

  ggsave(
    filename = file.path("plots", "umap", "08_umap_unintegrated_res0.4_split_samples.png"),
    plot = p3,
    width = 14,
    height = 8,
    dpi = 300
  )

  log_message("Stage 08_unintegrated_clustering_umap complete")
}

main()
