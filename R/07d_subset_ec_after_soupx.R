#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(yaml)
})

repo_root <- "/work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo"
setwd(repo_root)

decision_file <- "decisions/07_ec_retention_after_soupx.yml"
output_rds <- "results/objects/07_ec_enriched_after_soupx.rds"

dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

log_message <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", ts, paste(..., collapse = "")))
}

stop_if_missing <- function(path, label = NULL) {
  if (!file.exists(path) && !dir.exists(path)) {
    if (is.null(label)) {
      stop(sprintf("Required path does not exist: %s", path), call. = FALSE)
    } else {
      stop(sprintf("Required %s does not exist: %s", label, path), call. = FALSE)
    }
  }
}

log_message("Reading decision file")
stop_if_missing(decision_file, "decision file")
dec <- yaml::read_yaml(decision_file)

if (is.null(dec$source_object) || is.null(dec$retained_clusters) || is.null(dec$removed_clusters)) {
  stop("Decision file is missing one or more required fields: source_object, retained_clusters, removed_clusters.", call. = FALSE)
}

input_rds <- dec$source_object
retain_clusters <- as.character(unlist(dec$retained_clusters))
remove_clusters <- as.character(unlist(dec$removed_clusters))

log_message("Loading source object: ", input_rds)
stop_if_missing(input_rds, "source object")
obj <- readRDS(input_rds)

if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
  obj$seurat_clusters <- as.character(Idents(obj))
}

Idents(obj) <- "seurat_clusters"

present_clusters <- levels(Idents(obj))

missing_retain <- setdiff(retain_clusters, present_clusters)
if (length(missing_retain) > 0) {
  stop(
    sprintf(
      "Some retained clusters from decision file are absent from the object: %s",
      paste(missing_retain, collapse = ", ")
    ),
    call. = FALSE
  )
}

log_message("Retaining clusters: ", paste(retain_clusters, collapse = ", "))
log_message("Removing clusters: ", paste(remove_clusters, collapse = ", "))

obj_ec <- subset(
  x = obj,
  idents = retain_clusters
)

obj_ec$stage07d_source_cluster <- as.character(Idents(obj_ec))
Idents(obj_ec) <- "stage07d_source_cluster"

if (is.null(obj_ec@misc)) {
  obj_ec@misc <- list()
}

obj_ec@misc$stage07d_ec_subset <- list(
  decision_file = decision_file,
  source_object = input_rds,
  retained_clusters = retain_clusters,
  removed_clusters = remove_clusters,
  run_datetime = as.character(Sys.time())
)

saveRDS(obj_ec, output_rds)

log_message("Saved EC-enriched object: ", output_rds)
log_message("Number of retained cells: ", ncol(obj_ec))
log_message("Retained clusters present: ", paste(levels(Idents(obj_ec)), collapse = ", "))
log_message("Stage 07d completed successfully")
