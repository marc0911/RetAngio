suppressPackageStartupMessages({
  library(yaml)
})

timestamp_string <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

log_message <- function(...) {
  msg <- paste(..., collapse = " ")
  message(sprintf("%s | %s", timestamp_string(), msg))
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(normalizePath(path, winslash = "/", mustWork = FALSE))
}

assert_file_exists <- function(path, label = NULL) {
  if (is.null(path) || is.na(path) || !nzchar(path)) {
    stop(sprintf("%s is missing or empty.", ifelse(is.null(label), "File path", label)), call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("%s not found: %s", ifelse(is.null(label), "File", label), path), call. = FALSE)
  }
  invisible(path)
}

assert_dir_exists <- function(path, label = NULL) {
  if (is.null(path) || is.na(path) || !nzchar(path)) {
    stop(sprintf("%s is missing or empty.", ifelse(is.null(label), "Directory path", label)), call. = FALSE)
  }
  if (!dir.exists(path)) {
    stop(sprintf("%s not found: %s", ifelse(is.null(label), "Directory", label), path), call. = FALSE)
  }
  invisible(path)
}

read_config <- function(config_path = "config/config.yml") {
  assert_file_exists(config_path, "Config file")
  cfg <- yaml::read_yaml(config_path)

  required_top <- c("project", "input", "analysis")
  missing_top <- setdiff(required_top, names(cfg))
  if (length(missing_top) > 0) {
    stop(sprintf(
      "Config file is missing required top-level sections: %s",
      paste(missing_top, collapse = ", ")
    ), call. = FALSE)
  }

  cfg
}

save_rds_with_log <- function(object, path) {
  ensure_dir(dirname(path))
  saveRDS(object, path)
  log_message("Saved RDS:", path)
  invisible(path)
}

load_rds_with_log <- function(path) {
  assert_file_exists(path, "RDS file")
  log_message("Loading RDS:", path)
  readRDS(path)
}

read_decision_yaml <- function(path) {
  assert_file_exists(path, "Decision YAML")
  yaml::read_yaml(path)
}

require_decision_fields <- function(decision, fields, path) {
  missing_fields <- fields[!fields %in% names(decision)]
  if (length(missing_fields) > 0) {
    stop(sprintf(
      "Decision file %s is missing required fields: %s",
      path,
      paste(missing_fields, collapse = ", ")
    ), call. = FALSE)
  }

  empty_fields <- vapply(
    fields,
    function(x) {
      val <- decision[[x]]
      is.null(val) || (length(val) == 1 && is.na(val)) || (is.character(val) && !nzchar(val))
    },
    logical(1)
  )

  if (any(empty_fields)) {
    stop(sprintf(
      "Decision file %s has empty required fields: %s",
      path,
      paste(fields[empty_fields], collapse = ", ")
    ), call. = FALSE)
  }

  invisible(TRUE)
}

stop_for_missing_decision <- function(path, message = NULL) {
  base_msg <- sprintf("Required decision file is missing: %s", path)
  if (is.null(message)) {
    stop(base_msg, call. = FALSE)
  } else {
    stop(sprintf("%s\n%s", base_msg, message), call. = FALSE)
  }
}

create_stage_dirs <- function(cfg) {
  output_dir <- cfg$project$output_dir

  dirs <- c(
    "results",
    output_dir,
    file.path(output_dir, "objects"),
    file.path(output_dir, "qc"),
    file.path(output_dir, "doublets"),
    file.path(output_dir, "ambient"),
    file.path(output_dir, "unintegrated"),
    file.path(output_dir, "integration_compare"),
    file.path(output_dir, "annotation"),
    file.path(output_dir, "dge"),
    "plots",
    file.path("plots", "00_build"),
    "logs",
    "data",
    "data/metadata"
  )

  invisible(vapply(dirs, ensure_dir, character(1)))
}

normalize_optional_path <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(trimws(x))) {
    return(NA_character_)
  }
  trimws(x)
}

guess_sample_paths <- function(sample_id, extracted_dir) {
  base_dir <- file.path(extracted_dir, sample_id)

  # Search a few likely places, keeping this helper conservative.
  candidates_filtered <- c(
    file.path(base_dir, sample_id, "outs", "per_sample_outs", sample_id, "count", "sample_filtered_feature_bc_matrix.h5"),
    file.path(base_dir, sample_id, sample_id, "outs", "per_sample_outs", sample_id, "count", "sample_filtered_feature_bc_matrix.h5"),
    file.path(base_dir, "outs", "per_sample_outs", sample_id, "count", "sample_filtered_feature_bc_matrix.h5")
  )

  candidates_raw <- c(
    file.path(base_dir, sample_id, "outs", "multi", "count", "raw_feature_bc_matrix.h5"),
    file.path(base_dir, sample_id, sample_id, "outs", "multi", "count", "raw_feature_bc_matrix.h5"),
    file.path(base_dir, "outs", "multi", "count", "raw_feature_bc_matrix.h5")
  )

  candidates_mol <- c(
    file.path(base_dir, sample_id, "outs", "multi", "count", "raw_molecule_info.h5"),
    file.path(base_dir, sample_id, sample_id, "outs", "multi", "count", "raw_molecule_info.h5"),
    file.path(base_dir, "outs", "multi", "count", "raw_molecule_info.h5")
  )

  first_existing <- function(paths) {
    hit <- paths[file.exists(paths)]
    if (length(hit) == 0) NA_character_ else hit[[1]]
  }

  data.frame(
    sample_id = sample_id,
    filtered_h5 = first_existing(candidates_filtered),
    raw_h5 = first_existing(candidates_raw),
    molecule_info_h5 = first_existing(candidates_mol),
    stringsAsFactors = FALSE
  )
}

validate_sample_manifest <- function(samples_df) {
  required_cols <- c(
    "sample_id",
    "condition",
    "timepoint",
    "filtered_h5",
    "raw_h5",
    "molecule_info_h5"
  )

  missing_cols <- setdiff(required_cols, colnames(samples_df))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Sample manifest is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ), call. = FALSE)
  }

  if (anyDuplicated(samples_df$sample_id) > 0) {
    dup_ids <- unique(samples_df$sample_id[duplicated(samples_df$sample_id)])
    stop(sprintf(
      "Duplicated sample_id values found in sample manifest: %s",
      paste(dup_ids, collapse = ", ")
    ), call. = FALSE)
  }

  required_nonempty <- c("sample_id", "condition", "timepoint", "filtered_h5")
  for (col in required_nonempty) {
    vals <- samples_df[[col]]
    bad <- is.na(vals) | !nzchar(trimws(as.character(vals)))
    if (any(bad)) {
      stop(sprintf(
        "Column '%s' contains missing/empty values for sample(s): %s",
        col,
        paste(samples_df$sample_id[bad], collapse = ", ")
      ), call. = FALSE)
    }
  }

  samples_df$sample_id <- trimws(as.character(samples_df$sample_id))
  samples_df$condition <- trimws(as.character(samples_df$condition))
  samples_df$timepoint <- trimws(as.character(samples_df$timepoint))
  samples_df$filtered_h5 <- trimws(as.character(samples_df$filtered_h5))
  samples_df$raw_h5 <- vapply(samples_df$raw_h5, normalize_optional_path, character(1))
  samples_df$molecule_info_h5 <- vapply(samples_df$molecule_info_h5, normalize_optional_path, character(1))

  samples_df
}

read_samples_csv <- function(samples_csv) {
  assert_file_exists(samples_csv, "Sample manifest CSV")
  df <- read.csv(samples_csv, stringsAsFactors = FALSE, check.names = FALSE)
  validate_sample_manifest(df)
}
