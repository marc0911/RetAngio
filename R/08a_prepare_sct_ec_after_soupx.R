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
  input_rds <- "results/objects/07_ec_enriched_after_soupx.rds"
  output_merged <- "results/objects/08a_ec_after_soupx_sct_merged.rds"
  output_summary <- "results/sct/08a_ec_after_soupx_sct_per_sample_summary.csv"
  output_dir <- "results/sct"

  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 08a_prepare_sct_ec_after_soupx")

  future::plan("sequential")
  options(future.globals.maxSize = 8 * 1024^3)

  message_ts("future plan set to sequential")
  message_ts("future.globals.maxSize set to 8 GiB")

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  DefaultAssay(seu) <- "RNA"

  message_ts("Joining RNA layers before SCT preparation")
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])

  if (!"orig.ident" %in% colnames(seu@meta.data)) {
    stop("orig.ident not found in metadata.")
  }

  sample_ids <- unique(seu$orig.ident)
  sample_ids <- sort(as.character(sample_ids))
  message_ts("Samples detected:", paste(sample_ids, collapse = ", "))

  seu.list <- SplitObject(seu, split.by = "orig.ident")
  message_ts("Split object into", length(seu.list), "sample-level objects")

  sct_list <- list()
  summary_list <- list()

  for (sid in sample_ids) {
    message_ts("Running SCTransform for sample:", sid)

    obj <- seu.list[[sid]]
    DefaultAssay(obj) <- "RNA"

    obj <- SCTransform(
      obj,
      assay = "RNA",
      new.assay.name = "SCT",
      vst.flavor = "v2",
      verbose = TRUE
    )

    DefaultAssay(obj) <- "SCT"

    sample_out <- file.path(output_dir, paste0("08a_", sid, "_EC_after_SoupX_SCT.rds"))
    saveRDS(obj, sample_out)
    message_ts("Saved SCT object for sample:", sid)

    sct_list[[sid]] <- obj
    summary_list[[sid]] <- data.frame(
      sample_id = sid,
      n_cells = ncol(obj),
      n_genes_rna = nrow(obj[["RNA"]]),
      n_genes_sct = nrow(obj[["SCT"]]),
      stringsAsFactors = FALSE
    )
  }

  merged <- sct_list[[1]]
  if (length(sct_list) > 1) {
    merged <- merge(sct_list[[1]], y = sct_list[-1])
  }

  DefaultAssay(merged) <- "SCT"

  if (is.null(merged@misc)) {
    merged@misc <- list()
  }

  merged@misc$stage08a_prepare_sct <- list(
    source_object = input_rds,
    samples = sample_ids,
    assay = "RNA",
    new_assay = "SCT",
    vst_flavor = "v2",
    future_plan = "sequential",
    future_globals_maxSize = 8 * 1024^3,
    run_datetime = as.character(Sys.time())
  )

  saveRDS(merged, output_merged)
  message_ts("Saved merged SCT object:", output_merged)

  summary_df <- bind_rows(summary_list)
  write_csv(summary_df, output_summary)

  message_ts("Stage 08a_prepare_sct_ec_after_soupx complete")
}

main()
