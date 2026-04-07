#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(future)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_rds <- "results/objects/13_ec_sct_pca_unintegrated.rds"
  output_rds <- "results/objects/16_ec_rpca_integrated.rds"
  features_csv <- "results/integration/16_ec_rpca_features_used.csv"

  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/integration", recursive = TRUE, showWarnings = FALSE)

  future::plan("sequential")
  options(future.globals.maxSize = 8 * 1024^3)

  message_ts("Starting stage 16_integrate_ec_only_rpca")
  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  DefaultAssay(seu) <- "SCT"

  if (!"orig.ident" %in% colnames(seu@meta.data)) {
    stop("orig.ident not found in metadata")
  }
  if (!"pca" %in% names(seu@reductions)) {
    stop("pca reduction not found in object")
  }

  message_ts("Reconstructing sample-level SCT objects from PCA object")
  seu.list <- SplitObject(seu, split.by = "orig.ident")

  for (i in seq_along(seu.list)) {
    DefaultAssay(seu.list[[i]]) <- "SCT"
  }

  features <- SelectIntegrationFeatures(object.list = seu.list, nfeatures = 3000)
  write_csv(data.frame(feature = features), features_csv)
  message_ts("Selected integration features:", length(features))

  VariableFeatures(seu) <- features

  message_ts("Splitting SCT assay into layers by sample for IntegrateLayers")
  seu[["SCT"]] <- split(seu[["SCT"]], f = seu$orig.ident)

  seu <- IntegrateLayers(
    object = seu,
    method = RPCAIntegration,
    assay = "SCT",
    orig.reduction = "pca",
    new.reduction = "integrated.rpca",
    features = features,
    verbose = TRUE
  )

  saveRDS(seu, output_rds)
  message_ts("Saved integrated object:", output_rds)
  message_ts("Stage 16_integrate_ec_only_rpca complete")
}

main()
