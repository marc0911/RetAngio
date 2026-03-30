#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
})

source("R/00_utils.R")

args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1L]] else "config/config.yml"

log_message("Starting stage 0: build fresh merged Seurat object")
log_message("Using config file: ", config_path)

cfg <- read_config(config_path)
create_stage_dirs(cfg)

samples_csv_path <- cfg$input$samples_csv
samples <- read_samples_csv(samples_csv_path)

sample_objects <- list()
summary_rows <- vector("list", nrow(samples))

for (i in seq_len(nrow(samples))) {
  sample_id <- samples$sample_id[[i]]
  condition <- samples$condition[[i]]
  timepoint <- samples$timepoint[[i]]
  filtered_ref <- samples$filtered_h5[[i]]
  raw_h5 <- samples$raw_h5[[i]]
  molecule_info_h5 <- samples$molecule_info_h5[[i]]

  log_message("Processing sample: ", sample_id)

  # Validate secondary inputs listed in manifest.
  assert_file_exists(raw_h5, label = paste0(sample_id, " raw_h5"))
  assert_file_exists(molecule_info_h5, label = paste0(sample_id, " molecule_info_h5"))

  resolved <- resolve_filtered_input_path(filtered_ref)
  input_path <- resolved$path

  counts <- if (resolved$type == "matrix_dir") {
    log_message("Reading 10X matrix directory for ", sample_id, ": ", input_path)
    Seurat::Read10X(data.dir = input_path)
  } else {
    log_message("Reading 10X H5 matrix for ", sample_id, ": ", input_path)
    Seurat::Read10X_h5(filename = input_path)
  }

  # If Read10X* returns a list (e.g., multi-assay), prefer "Gene Expression" when available.
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      assay_names <- names(counts)
      counts <- counts[[1L]]
      log_message(
        "Read10X returned multiple matrices for ", sample_id,
        "; using first entry: ", assay_names[[1L]]
      )
    }
  }

  so <- Seurat::CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = cfg$preprocessing$min_cells_per_gene,
    min.features = cfg$preprocessing$min_features_create_object
  )

  n_cells <- ncol(so)
  if (n_cells <= 0L) {
    stop(sprintf("No cells loaded for sample '%s'. Aborting.", sample_id), call. = FALSE)
  }

  so$sample_id <- sample_id
  so$condition <- condition
  so$timepoint <- timepoint

  # Keep file-level provenance at object level (cleaner than repeating per-cell paths).
  so@misc$sample_manifest <- list(
    sample_id = sample_id,
    condition = condition,
    timepoint = timepoint,
    filtered_input = filtered_ref,
    filtered_input_resolved = input_path,
    raw_h5 = raw_h5,
    molecule_info_h5 = molecule_info_h5
  )

  # Explicitly align orig.ident with sample_id for downstream grouping consistency.
  so$orig.ident <- sample_id

  sample_objects[[sample_id]] <- so
  summary_rows[[i]] <- data.frame(
    sample_id = sample_id,
    n_cells = n_cells,
    condition = condition,
    timepoint = timepoint,
    stringsAsFactors = FALSE
  )

  log_message("Loaded sample ", sample_id, " with ", n_cells, " cells")
}

if (length(sample_objects) == 0L) {
  stop("No sample objects were created. Check config/samples.csv.", call. = FALSE)
}

log_message("Merging ", length(sample_objects), " sample objects")
sample_ids <- names(sample_objects)
merged_object <- if (length(sample_objects) == 1L) {
  sample_objects[[1L]]
} else {
  Seurat::merge(
    x = sample_objects[[1L]],
    y = sample_objects[2:length(sample_objects)],
    add.cell.ids = sample_ids,
    project = cfg$project$name
  )
}

merged_object@misc$sample_manifest <- samples

summary_df <- do.call(rbind, summary_rows)

save_rds_with_log(merged_object, "results/objects/00_raw_merged.rds")
utils::write.csv(samples, "data/metadata/samples_resolved.csv", row.names = FALSE)
log_message("Saved resolved sample manifest: data/metadata/samples_resolved.csv")
utils::write.csv(summary_df, "results/qc/00_build_summary.csv", row.names = FALSE)
log_message("Saved build summary: results/qc/00_build_summary.csv")

log_message("Stage 0 completed successfully")
