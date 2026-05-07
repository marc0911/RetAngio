#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

repo_root <- "/work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo"
setwd(repo_root)

input_rds <- "results/objects/17_ec_integrated_clustered_umap.rds"
out_txt <- "logs/19b_inspect_plotting_metadata.txt"

dir.create("logs", recursive = TRUE, showWarnings = FALSE)

log_message <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", ts, paste(..., collapse = "")))
}

capture_block <- function(..., file, append = TRUE) {
  cat(..., file = file, append = append, sep = "")
}

log_message("Loading object: ", input_rds)
obj <- readRDS(input_rds)
md <- obj@meta.data

if (file.exists(out_txt)) {
  file.remove(out_txt)
}

capture_block("=== METADATA COLUMN NAMES ===\n", file = out_txt, append = TRUE)
capture_block(paste(colnames(md), collapse = "\n"), "\n\n", file = out_txt, append = TRUE)

candidate_cols <- c(
  "sample_id", "condition", "group", "timepoint", "orig.ident",
  "sample", "dataset", "batch", "age", "stage"
)

present_cols <- candidate_cols[candidate_cols %in% colnames(md)]

capture_block("=== CANDIDATE PLOTTING COLUMNS PRESENT ===\n", file = out_txt, append = TRUE)
if (length(present_cols) == 0) {
  capture_block("None\n\n", file = out_txt, append = TRUE)
} else {
  capture_block(paste(present_cols, collapse = "\n"), "\n\n", file = out_txt, append = TRUE)
}

for (cc in present_cols) {
  vals <- unique(as.character(md[[cc]]))
  vals <- vals[order(vals)]
  vals <- vals[!is.na(vals)]

  capture_block(sprintf("=== UNIQUE VALUES: %s ===\n", cc), file = out_txt, append = TRUE)
  if (length(vals) == 0) {
    capture_block("No non-NA values\n\n", file = out_txt, append = TRUE)
  } else {
    capture_block(paste(vals, collapse = "\n"), "\n\n", file = out_txt, append = TRUE)
  }
}

capture_block("=== FIRST 20 ROWS OF KEY COLUMNS ===\n", file = out_txt, append = TRUE)
if (length(present_cols) > 0) {
  preview_df <- md[, present_cols, drop = FALSE]
  preview_df <- head(preview_df, 20)
  capture.output(print(preview_df), file = out_txt, append = TRUE)
  capture_block("\n", file = out_txt, append = TRUE)
}

if (!is.null(obj@misc)) {
  capture_block("=== MISC NAMES ===\n", file = out_txt, append = TRUE)
  capture_block(paste(names(obj@misc), collapse = "\n"), "\n\n", file = out_txt, append = TRUE)

  if (!is.null(obj@misc$sample_manifest) && is.data.frame(obj@misc$sample_manifest)) {
    capture_block("=== MISC$sample_manifest COLUMNS ===\n", file = out_txt, append = TRUE)
    capture_block(paste(colnames(obj@misc$sample_manifest), collapse = "\n"), "\n\n", file = out_txt, append = TRUE)

    capture_block("=== FIRST 20 ROWS OF MISC$sample_manifest ===\n", file = out_txt, append = TRUE)
    capture.output(print(utils::head(obj@misc$sample_manifest, 20)), file = out_txt, append = TRUE)
    capture_block("\n", file = out_txt, append = TRUE)
  }
}

log_message("Metadata inspection written to: ", out_txt)
