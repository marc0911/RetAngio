# Utility helpers for RetAngio stage scripts.
# Keep this file lightweight and dependency-minimal.

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(normalizePath(path, winslash = "/", mustWork = FALSE))
}

log_message <- function(...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", timestamp, paste(..., collapse = "")))
}

assert_file_exists <- function(path, label = NULL) {
  if (!file.exists(path)) {
    if (is.null(label)) {
      stop(sprintf("Required file does not exist: '%s'", path), call. = FALSE)
    }
    stop(sprintf("Required file (%s) does not exist: '%s'", label, path), call. = FALSE)
  }
  invisible(TRUE)
}

assert_dir_exists <- function(path, label = NULL) {
  if (!dir.exists(path)) {
    if (is.null(label)) {
      stop(sprintf("Required directory does not exist: '%s'", path), call. = FALSE)
    }
    stop(sprintf("Required directory (%s) does not exist: '%s'", label, path), call. = FALSE)
  }
  invisible(TRUE)
}

read_config <- function(config_path) {
  assert_file_exists(config_path, label = "config")
  cfg <- yaml::read_yaml(config_path)
  if (!is.list(cfg)) {
    stop(sprintf("Config file '%s' did not parse as a YAML mapping/list.", config_path), call. = FALSE)
  }
  cfg
}

save_rds_with_log <- function(object, path) {
  ensure_dir(dirname(path))
  saveRDS(object, path)
  log_message("Saved RDS: ", path)
}

load_rds_with_log <- function(path) {
  assert_file_exists(path, label = "RDS")
  log_message("Loading RDS: ", path)
  readRDS(path)
}

read_samples_csv <- function(samples_csv) {
  assert_file_exists(samples_csv, label = "samples CSV")

  samples <- utils::read.csv(
    file = samples_csv,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  required_cols <- c(
    "sample_id", "condition", "timepoint",
    "filtered_h5", "raw_h5", "molecule_info_h5"
  )
  missing_cols <- setdiff(required_cols, colnames(samples))
  if (length(missing_cols) > 0L) {
    stop(
      sprintf(
        "Sample manifest is missing required column(s): %s",
        paste(missing_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (anyDuplicated(samples$sample_id) > 0L) {
    dup_ids <- unique(samples$sample_id[duplicated(samples$sample_id)])
    stop(
      sprintf("Duplicate sample_id values found: %s", paste(dup_ids, collapse = ", ")),
      call. = FALSE
    )
  }

  required_nonempty_stage0 <- c("sample_id", "condition", "timepoint", "filtered_h5")
  for (col in required_nonempty_stage0) {
    is_missing <- is.na(samples[[col]]) | trimws(samples[[col]]) == ""
    if (any(is_missing)) {
      bad_rows <- which(is_missing)
      stop(
        sprintf(
          "Stage-0 required column '%s' has missing value(s) at row(s): %s",
          col,
          paste(bad_rows, collapse = ", ")
        ),
        call. = FALSE
      )
    }
  }

  # raw_h5 and molecule_info_h5 are optional for stage 0.
  # Keep as character, normalize blanks to NA for cleaner downstream handling.
  optional_cols <- c("raw_h5", "molecule_info_h5")
  for (col in optional_cols) {
    vals <- trimws(samples[[col]])
    vals[is.na(vals) | vals == ""] <- NA_character_
    samples[[col]] <- vals
  }

  samples
}

create_stage_dirs <- function(cfg) {
  dirs <- c(
    cfg$project$output_dir,
    "results/objects",
    "results/qc",
    "results/doublets",
    "results/ambient",
    "results/unintegrated",
    "results/integration_compare",
    "results/annotation",
    "results/dge",
    "plots",
    "logs",
    "data/metadata"
  )
  invisible(lapply(dirs, ensure_dir))
}

read_decision_yaml <- function(path) {
  assert_file_exists(path, label = "decision file")
  decision <- yaml::read_yaml(path)
  if (!is.list(decision)) {
    stop(sprintf("Decision file '%s' is not a valid YAML mapping/list.", path), call. = FALSE)
  }
  decision
}

require_decision_fields <- function(decision, fields, path) {
  missing <- setdiff(fields, names(decision))
  if (length(missing) > 0L) {
    stop(
      sprintf(
        "Decision file '%s' missing field(s): %s",
        path,
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

stop_for_missing_decision <- function(path, message = NULL) {
  if (!file.exists(path)) {
    if (is.null(message)) {
      message <- sprintf("Missing required decision file: '%s'", path)
    }
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

guess_sample_paths <- function(sample_id, extracted_dir = "data/extracted") {
  # Convenience helper only. The manifest remains source of truth.
  sample_dir <- file.path(extracted_dir, sample_id)
  if (!dir.exists(sample_dir)) {
    return(list(filtered = NA_character_, raw = NA_character_, molecule = NA_character_))
  }

  filtered_candidates <- c(
    file.path(sample_dir, "sample_filtered_feature_bc_matrix.h5"),
    file.path(sample_dir, "filtered_feature_bc_matrix.h5"),
    file.path(sample_dir, "sample_filtered_feature_bc_matrix")
  )
  raw_candidates <- c(
    file.path(sample_dir, "raw_feature_bc_matrix.h5"),
    file.path(sample_dir, "raw_feature_bc_matrix")
  )
  molecule_candidates <- c(
    file.path(sample_dir, "raw_molecule_info.h5"),
    file.path(sample_dir, "molecule_info.h5")
  )

  pick_existing <- function(paths) {
    idx <- which(file.exists(paths) | dir.exists(paths))
    if (length(idx) == 0L) {
      return(NA_character_)
    }
    paths[[idx[[1L]]]]
  }

  list(
    filtered = pick_existing(filtered_candidates),
    raw = pick_existing(raw_candidates),
    molecule = pick_existing(molecule_candidates)
  )
}

resolve_filtered_input_path <- function(filtered_path) {
  # Primary expectation: filtered matrix path from manifest.
  if (file.exists(filtered_path)) {
    return(list(path = filtered_path, type = "h5_or_file"))
  }

  if (dir.exists(filtered_path)) {
    return(list(path = filtered_path, type = "matrix_dir"))
  }

  # Conservative fallback if a sample directory was supplied instead of matrix path.
  dir_candidate <- file.path(filtered_path, "sample_filtered_feature_bc_matrix")
  if (dir.exists(dir_candidate)) {
    return(list(path = dir_candidate, type = "matrix_dir"))
  }

  h5_candidates <- c(
    filtered_path,
    file.path(filtered_path, "sample_filtered_feature_bc_matrix.h5"),
    file.path(filtered_path, "filtered_feature_bc_matrix.h5")
  )

  for (candidate in h5_candidates) {
    if (file.exists(candidate)) {
      return(list(path = candidate, type = "h5_or_file"))
    }
  }

  stop(
    sprintf(
      paste(
        "Unable to resolve filtered matrix input from manifest value '%s'.",
        "Provide a valid filtered H5 path or filtered matrix directory."
      ),
      filtered_path
    ),
    call. = FALSE
  )
}
