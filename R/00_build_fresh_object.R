args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1) args[[1]] else "config/config.yml"

source("R/00_utils.R")

suppressPackageStartupMessages({
  library(Seurat)
})

read_filtered_matrix <- function(filtered_h5_path) {
  # Preferred route: direct H5 import from Cell Ranger.
  if (!is.na(filtered_h5_path) && file.exists(filtered_h5_path)) {
    mat <- Seurat::Read10X_h5(filtered_h5_path)
    if (is.list(mat)) {
      # Keep Gene Expression if present; otherwise keep first entry.
      if ("Gene Expression" %in% names(mat)) {
        mat <- mat[["Gene Expression"]]
      } else {
        mat <- mat[[1]]
      }
    }
    return(mat)
  }

  stop(sprintf("Filtered matrix file not found: %s", filtered_h5_path), call. = FALSE)
}

main <- function() {
  log_message("Starting stage 00_build_fresh_object")
  log_message("Using config:", config_path)

  cfg <- read_config(config_path)
  create_stage_dirs(cfg)

  set.seed(cfg$project$seed)

  samples_df <- read_samples_csv(cfg$input$samples_csv)

  log_message("Loaded sample manifest with", nrow(samples_df), "samples")

  # Stage 0 requires filtered matrix input.
  missing_filtered <- !file.exists(samples_df$filtered_h5)
  if (any(missing_filtered)) {
    stop(sprintf(
      "Missing filtered_h5 for sample(s): %s",
      paste(samples_df$sample_id[missing_filtered], collapse = ", ")
    ), call. = FALSE)
  }

  # raw_h5 and molecule_info_h5 are optional at stage 0; keep them for later stages.
  missing_raw <- is.na(samples_df$raw_h5) | !file.exists(samples_df$raw_h5)
  missing_mol <- is.na(samples_df$molecule_info_h5) | !file.exists(samples_df$molecule_info_h5)

  if (any(missing_raw)) {
    log_message(
      "Warning: raw_h5 missing for sample(s):",
      paste(samples_df$sample_id[missing_raw], collapse = ", "),
      "- continuing because stage 0 only requires filtered matrices."
    )
    samples_df$raw_h5[missing_raw] <- NA_character_
  }

  if (any(missing_mol)) {
    log_message(
      "Warning: molecule_info_h5 missing for sample(s):",
      paste(samples_df$sample_id[missing_mol], collapse = ", "),
      "- continuing because stage 0 only requires filtered matrices."
    )
    samples_df$molecule_info_h5[missing_mol] <- NA_character_
  }

  seu_list <- vector("list", length = nrow(samples_df))
  names(seu_list) <- samples_df$sample_id
  summary_rows <- vector("list", length = nrow(samples_df))

  for (i in seq_len(nrow(samples_df))) {
    sample_row <- samples_df[i, , drop = FALSE]
    sample_id <- sample_row$sample_id

    log_message("Reading sample:", sample_id)
    counts <- read_filtered_matrix(sample_row$filtered_h5)

    if (ncol(counts) == 0) {
      stop(sprintf("No cells loaded for sample %s", sample_id), call. = FALSE)
    }

    seu <- CreateSeuratObject(
      counts = counts,
      project = sample_id,
      min.cells = cfg$preprocessing$min_cells_per_gene,
      min.features = cfg$preprocessing$min_features_create_object
    )

    seu$sample_id <- sample_row$sample_id
    seu$condition <- sample_row$condition
    seu$timepoint <- sample_row$timepoint
    seu$orig.ident <- sample_row$sample_id

    # Cleaner than repeating file paths in per-cell metadata.
    seu@misc$sample_manifest <- list(
      sample_id = sample_row$sample_id,
      condition = sample_row$condition,
      timepoint = sample_row$timepoint,
      filtered_h5 = sample_row$filtered_h5,
      raw_h5 = sample_row$raw_h5,
      molecule_info_h5 = sample_row$molecule_info_h5
    )

    seu_list[[i]] <- seu

    summary_rows[[i]] <- data.frame(
      sample_id = sample_row$sample_id,
      n_cells = ncol(seu),
      condition = sample_row$condition,
      timepoint = sample_row$timepoint,
      stringsAsFactors = FALSE
    )

    log_message("Loaded", ncol(seu), "cells for", sample_id)
  }

  build_summary <- do.call(rbind, summary_rows)

  if (length(seu_list) == 1) {
    merged_seu <- seu_list[[1]]
  } else {
    merged_seu <- merge(
      x = seu_list[[1]],
      y = seu_list[-1],
      add.cell.ids = names(seu_list),
      project = cfg$project$name
    )
  }

  # Keep the full resolved sample manifest at object level for downstream stages.
  merged_seu@misc$sample_manifest <- samples_df
  merged_seu@misc$build_stage <- "00_raw_merged"
  merged_seu@misc$config_path <- config_path

  save_rds_with_log(
    merged_seu,
    file.path(cfg$project$output_dir, "objects", "00_raw_merged.rds")
  )

  write.csv(
    samples_df,
    file = file.path("data", "metadata", "samples_resolved.csv"),
    row.names = FALSE,
    quote = TRUE
  )
  log_message("Saved resolved sample manifest: data/metadata/samples_resolved.csv")

  write.csv(
    build_summary,
    file = file.path(cfg$project$output_dir, "qc", "00_build_summary.csv"),
    row.names = FALSE,
    quote = TRUE
  )
  log_message("Saved build summary:", file.path(cfg$project$output_dir, "qc", "00_build_summary.csv"))

  log_message("Stage 00_build_fresh_object complete")
}

main()
