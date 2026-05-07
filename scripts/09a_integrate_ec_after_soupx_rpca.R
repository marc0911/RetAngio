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
  input_rds <- "results/objects/08b_ec_after_soupx_sct_pca_unintegrated.rds"
  output_rds <- "results/objects/09a_ec_after_soupx_rpca_integrated.rds"
  features_csv <- "results/integration/09a_ec_after_soupx_rpca_features_used.csv"

  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/integration", recursive = TRUE, showWarnings = FALSE)

  future::plan("sequential")
  options(future.globals.maxSize = 8 * 1024^3)

  message_ts("Starting stage 09a_integrate_ec_after_soupx_rpca")
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

  if (is.null(seu@misc)) {
    seu@misc <- list()
  }

  seu@misc$stage09a_rpca_integration <- list(
    source_object = input_rds,
    method = "RPCAIntegration",
    assay = "SCT",
    orig_reduction = "pca",
    new_reduction = "integrated.rpca",
    nfeatures = length(features),
    run_datetime = as.character(Sys.time())
  )

  saveRDS(seu, output_rds)
  message_ts("Saved integrated object:", output_rds)
  message_ts("Stage 09a_integrate_ec_after_soupx_rpca complete")
}

main()
