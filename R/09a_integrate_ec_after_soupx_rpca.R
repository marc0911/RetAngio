#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(dplyr)
  library(readr)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_merged <- "results/objects/08a_ec_after_soupx_sct_merged.rds"
  input_sct_dir <- "results/sct"
  output_rds <- "results/objects/09a_ec_after_soupx_rpca_integrated.rds"
  output_summary <- "results/integration/09a_ec_after_soupx_rpca_summary.csv"

  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/integration", recursive = TRUE, showWarnings = FALSE)

  future::plan("sequential")
  options(future.globals.maxSize = 8 * 1024^3)

  message_ts("Starting stage 09a_integrate_ec_after_soupx_rpca")
  message_ts("future plan set to sequential")
  message_ts("future.globals.maxSize set to 8 GiB")

  merged <- readRDS(input_merged)
  message_ts("Loaded merged object:", input_merged)

  DefaultAssay(merged) <- "SCT"

  if (!"orig.ident" %in% colnames(merged@meta.data)) {
    stop("orig.ident not found in metadata.")
  }

  sample_ids <- sort(unique(as.character(merged$orig.ident)))
  message_ts("Samples detected:", paste(sample_ids, collapse = ", "))

  sct_paths <- file.path(input_sct_dir, paste0("08a_", sample_ids, "_EC_after_SoupX_SCT.rds"))
  missing_paths <- sct_paths[!file.exists(sct_paths)]
  if (length(missing_paths) > 0) {
    stop(paste("Missing per-sample SCT objects:", paste(missing_paths, collapse = ", ")))
  }

  message_ts("Loading per-sample SCT objects from results/sct")
  seu.list <- lapply(sct_paths, readRDS)
  names(seu.list) <- sample_ids

  for (sid in names(seu.list)) {
    DefaultAssay(seu.list[[sid]]) <- "SCT"
  }

  message_ts("Selecting integration features from per-sample SCT objects")
  features <- SelectIntegrationFeatures(
    object.list = seu.list,
    nfeatures = 3000
  )

  # Conservative filter: keep only features present in every object
  common_features <- Reduce(intersect, lapply(seu.list, rownames))
  features <- features[features %in% common_features]
  message_ts("Integration features after common-feature filter:", length(features))

  if (length(features) < 1000) {
    stop("Too few shared integration features remain after filtering.", call. = FALSE)
  }

  message_ts("Running PrepSCTIntegration")
  seu.list <- PrepSCTIntegration(
    object.list = seu.list,
    anchor.features = features,
    verbose = TRUE
  )
  message_ts("PrepSCTIntegration complete")

  for (sid in names(seu.list)) {
    message_ts("Running PCA for sample:", sid)
    DefaultAssay(seu.list[[sid]]) <- "SCT"
    VariableFeatures(seu.list[[sid]]) <- features
    seu.list[[sid]] <- RunPCA(
      object = seu.list[[sid]],
      assay = "SCT",
      features = features,
      npcs = 30,
      verbose = TRUE
    )
  }

  message_ts("Finding RPCA anchors")
  anchors <- FindIntegrationAnchors(
    object.list = seu.list,
    normalization.method = "SCT",
    anchor.features = features,
    reduction = "rpca",
    dims = 1:30,
    verbose = TRUE
  )

  message_ts("Integrating data")
  integrated <- IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT",
    dims = 1:30,
    verbose = TRUE
  )

  DefaultAssay(integrated) <- "integrated"

  if (is.null(integrated@misc)) {
    integrated@misc <- list()
  }

  integrated@misc$stage09a_rpca_integration <- list(
    source_object = input_merged,
    source_sct_dir = input_sct_dir,
    samples = sample_ids,
    normalization_method = "SCT",
    reduction = "rpca",
    dims_used = 1:30,
    nfeatures = length(features),
    run_datetime = as.character(Sys.time())
  )

  saveRDS(integrated, output_rds)

  summary_df <- data.frame(
    sample_id = sample_ids,
    n_cells = vapply(seu.list, ncol, numeric(1)),
    stringsAsFactors = FALSE
  )
  write_csv(summary_df, output_summary)

  message_ts("Saved RPCA-integrated object:", output_rds)
  message_ts("Saved RPCA summary:", output_summary)
  message_ts("Stage 09a_integrate_ec_after_soupx_rpca complete")
}

main()
