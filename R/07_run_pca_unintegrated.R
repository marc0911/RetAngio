source("R/00_utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tibble)
})

main <- function() {
  log_message("Starting stage 07_run_pca_unintegrated")

  cfg <- read_config()
  outdir <- cfg$project$output_dir
  sample_col <- cfg$analysis$sample_column

  ensure_dir(outdir)
  ensure_dir(file.path(outdir, "objects"))
  ensure_dir(file.path(outdir, "pca"))
  ensure_dir("plots")
  ensure_dir(file.path("plots", "pca"))

  in_rds <- file.path(outdir, "objects", "06_sct_merged.rds")
  if (!file.exists(in_rds)) {
    stop("Missing input object: ", in_rds)
  }

  seu <- readRDS(in_rds)
  log_message("Loaded object: ", in_rds)

  if (!"SCT" %in% names(seu@assays)) {
    stop("SCT assay not found in object.")
  }
  if (!sample_col %in% colnames(seu@meta.data)) {
    stop("Sample column not found in metadata: ", sample_col)
  }

  DefaultAssay(seu) <- "SCT"

  # Recreate sample-level list to define robust shared variable features
  seu_list <- SplitObject(seu, split.by = sample_col)
  log_message("Reconstructed ", length(seu_list), " sample-level objects from merged SCT object")

  nfeatures_use <- 3000
  features_use <- SelectIntegrationFeatures(
    object.list = seu_list,
    nfeatures = nfeatures_use
  )

  VariableFeatures(seu) <- features_use
  log_message("Selected ", length(features_use), " shared variable features for PCA")

  npcs <- cfg$unintegrated$npcs
  log_message("Running PCA with npcs = ", npcs)

  seu <- RunPCA(
    object = seu,
    assay = "SCT",
    features = features_use,
    npcs = npcs,
    verbose = TRUE
  )

  out_pca <- file.path(outdir, "objects", "07_sct_pca_unintegrated.rds")
  saveRDS(seu, out_pca)
  log_message("Saved PCA object: ", out_pca)

  write_csv(
    tibble(feature = features_use),
    file.path(outdir, "pca", "07_pca_features_used.csv")
  )

  sdev <- seu[["pca"]]@stdev
  pc_tbl <- tibble(
    PC = seq_along(sdev),
    stdev = sdev,
    variance = sdev^2,
    pct_variance = 100 * (variance / sum(variance)),
    cum_pct_variance = cumsum(pct_variance)
  )

  write_csv(pc_tbl, file.path(outdir, "pca", "07_pca_variance_summary.csv"))
  log_message("Saved PCA variance summary")

  p_elbow <- ElbowPlot(seu, ndims = npcs)

  ggsave(
    filename = file.path("plots", "pca", "07_elbowplot_unintegrated.png"),
    plot = p_elbow,
    width = 7,
    height = 5,
    dpi = 300
  )

  p_var <- ggplot(pc_tbl, aes(x = PC, y = pct_variance)) +
    geom_line() +
    geom_point() +
    theme_bw()

  ggsave(
    filename = file.path("plots", "pca", "07_pca_pct_variance.png"),
    plot = p_var,
    width = 7,
    height = 5,
    dpi = 300
  )

  p_cum <- ggplot(pc_tbl, aes(x = PC, y = cum_pct_variance)) +
    geom_line() +
    geom_point() +
    theme_bw()

  ggsave(
    filename = file.path("plots", "pca", "07_pca_cumulative_variance.png"),
    plot = p_cum,
    width = 7,
    height = 5,
    dpi = 300
  )

  log_message("Stage 07_run_pca_unintegrated complete")
}

main()
