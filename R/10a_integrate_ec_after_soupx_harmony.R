#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(future)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_rds <- "results/objects/08b_ec_after_soupx_sct_pca_unintegrated.rds"
  output_rds <- "results/objects/10a_ec_after_soupx_harmony_integrated.rds"

  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)

  future::plan("sequential")
  options(future.globals.maxSize = 8 * 1024^3)

  message_ts("Starting stage 10a_integrate_ec_after_soupx_harmony")
  message_ts("future plan set to sequential")
  message_ts("future.globals.maxSize set to 8 GiB")

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  DefaultAssay(seu) <- "SCT"

  if (!"orig.ident" %in% colnames(seu@meta.data)) {
    stop("orig.ident not found in metadata")
  }
  if (!"pca" %in% names(seu@reductions)) {
    stop("pca reduction not found in object")
  }

  message_ts("Running Harmony on PCA reduction using orig.ident")
  seu <- RunHarmony(
    object = seu,
    group.by.vars = "orig.ident",
    reduction = "pca",
    assay.use = "SCT",
    reduction.save = "harmony",
    dims.use = 1:30,
    verbose = TRUE
  )

  if (!"harmony" %in% names(seu@reductions)) {
    stop("Harmony reduction was not created")
  }

  if (is.null(seu@misc)) {
    seu@misc <- list()
  }

  seu@misc$stage10a_harmony_integration <- list(
    source_object = input_rds,
    method = "RunHarmony",
    assay = "SCT",
    input_reduction = "pca",
    output_reduction = "harmony",
    group_by = "orig.ident",
    dims_used = 1:30,
    run_datetime = as.character(Sys.time())
  )

  saveRDS(seu, output_rds)

  message_ts("Saved Harmony-integrated object:", output_rds)
  message_ts("Stage 10a_integrate_ec_after_soupx_harmony complete")
}

main()
