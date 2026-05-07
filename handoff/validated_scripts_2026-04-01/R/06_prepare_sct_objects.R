source("R/00_utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
  library(future)
})

main <- function() {
  log_message("Starting stage 06_prepare_sct_objects")

  cfg <- read_config()
  outdir <- cfg$project$output_dir
  sample_col <- cfg$analysis$sample_column

  ensure_dir(outdir)
  ensure_dir(file.path(outdir, "objects"))
  ensure_dir(file.path(outdir, "sct"))

  in_rds <- file.path(outdir, "objects", "04_qc_doublet_filtered.rds")
  if (!file.exists(in_rds)) {
    stop("Missing input object: ", in_rds)
  }

  seu <- readRDS(in_rds)
  log_message("Loaded object: ", in_rds)

  if (!sample_col %in% colnames(seu@meta.data)) {
    stop("Sample column not found in metadata: ", sample_col)
  }

  DefaultAssay(seu) <- "RNA"

  if ("RNA" %in% names(seu@assays)) {
    log_message("Joining RNA layers before SCT preparation")
    seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  }

  sample_ids <- unique(as.character(seu@meta.data[[sample_col]]))
  sample_ids <- sort(sample_ids)
  log_message("Samples detected: ", paste(sample_ids, collapse = ", "))

  seu_list <- SplitObject(seu, split.by = sample_col)
  log_message("Split object into ", length(seu_list), " sample-level objects")

  # Robust setting for SCTransform on cluster
  future::plan("sequential")
  options(future.globals.maxSize = 8 * 1024^3)

  summary_list <- vector("list", length(seu_list))
  names(summary_list) <- names(seu_list)

  for (sid in names(seu_list)) {
    obj <- seu_list[[sid]]
    DefaultAssay(obj) <- "RNA"

    log_message("Running SCTransform for sample: ", sid)

    obj <- SCTransform(
      obj,
      assay = "RNA",
      new.assay.name = "SCT",
      vst.flavor = "v2",
      verbose = TRUE
    )

    saveRDS(obj, file.path(outdir, "sct", paste0("06_", sid, "_SCT.rds")))
    log_message("Saved SCT object for sample: ", sid)

    seu_list[[sid]] <- obj

    summary_list[[sid]] <- tibble(
      sample_id = sid,
      n_cells = ncol(obj),
      n_genes_rna = nrow(obj[["RNA"]]),
      n_genes_sct = nrow(obj[["SCT"]])
    )
  }

  summary_tbl <- bind_rows(summary_list)
  write_csv(summary_tbl, file.path(outdir, "sct", "06_sct_per_sample_summary.csv"))

  merged_sct <- merge(
    x = seu_list[[1]],
    y = seu_list[-1]
  )

  saveRDS(merged_sct, file.path(outdir, "objects", "06_sct_merged.rds"))
  log_message("Saved merged SCT object: ", file.path(outdir, "objects", "06_sct_merged.rds"))

  log_message("Stage 06_prepare_sct_objects complete")
}

main()
