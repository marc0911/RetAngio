#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
  library(Matrix)
})

repo_root <- "/work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo"
setwd(repo_root)

input_rds <- "results/objects/05a_qc_doublet_precluster_for_soupx.rds"
output_rds <- "results/objects/05_soupx_corrected_merged.rds"
output_csv <- "results/ambient/05_soupx_contamination_summary.csv"

dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
dir.create("results/ambient", recursive = TRUE, showWarnings = FALSE)
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

read_10x_flexible <- function(path_value) {
  if (is.na(path_value) || trimws(path_value) == "") {
    stop("Encountered empty matrix path in sample manifest.", call. = FALSE)
  }

  path_value <- trimws(path_value)

  if (file.exists(path_value)) {
    if (grepl("\\.h5$", path_value, ignore.case = TRUE)) {
      x <- Read10X_h5(path_value)
    } else {
      stop(sprintf("Unsupported file path: %s", path_value), call. = FALSE)
    }
  } else if (dir.exists(path_value)) {
    x <- Read10X(data.dir = path_value)
  } else {
    stop(sprintf("Matrix path does not exist: %s", path_value), call. = FALSE)
  }

  if (is.list(x)) {
    if ("Gene Expression" %in% names(x)) {
      x <- x[["Gene Expression"]]
    } else {
      x <- x[[1]]
    }
  }

  if (!inherits(x, "dgCMatrix")) {
    x <- as(x, "dgCMatrix")
  }

  x
}

extract_sample_manifest <- function(obj) {
  if (!is.null(obj@misc$sample_manifest)) {
    manifest <- obj@misc$sample_manifest
    if (is.data.frame(manifest) && nrow(manifest) > 0) {
      return(manifest)
    }
  }

  cfg <- "config/samples.csv"
  stop_if_missing(cfg, "sample manifest")
  read.csv(cfg, stringsAsFactors = FALSE, check.names = FALSE)
}

pick_manifest_col <- function(df, candidates, label) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) == 0) {
    stop(
      sprintf(
        "Could not find a '%s' column in sample manifest. Tried: %s",
        label,
        paste(candidates, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  hit[[1]]
}

extract_raw_barcodes <- function(cell_names, sample_id) {
  prefix <- paste0(sample_id, "_")
  has_prefix <- startsWith(cell_names, prefix)
  if (!all(has_prefix)) {
    stop(
      sprintf("Some cells for sample %s do not start with expected prefix '%s'", sample_id, prefix),
      call. = FALSE
    )
  }
  sub(prefix, "", cell_names, fixed = TRUE)
}

log_message("Loading preclustered post-doublet object")
stop_if_missing(input_rds, "preclustered post-doublet object")
obj <- readRDS(input_rds)

if (!"sample_id" %in% colnames(obj@meta.data)) {
  stop("The input object must contain a 'sample_id' column in metadata.", call. = FALSE)
}

if (!"pre_soupx_cluster" %in% colnames(obj@meta.data)) {
  stop("The input object must contain a 'pre_soupx_cluster' column in metadata.", call. = FALSE)
}

manifest <- extract_sample_manifest(obj)

sample_col <- pick_manifest_col(manifest, c("sample_id", "sample", "SampleID"), "sample identifier")
filtered_col <- pick_manifest_col(
  manifest,
  c("filtered_h5", "filtered_matrix_path", "filtered_path"),
  "filtered matrix path"
)
raw_col <- pick_manifest_col(
  manifest,
  c("raw_h5", "raw_matrix_path", "raw_path"),
  "raw matrix path"
)

manifest <- manifest[, c(sample_col, filtered_col, raw_col)]
colnames(manifest) <- c("sample_id", "filtered_path", "raw_path")

sample_ids <- unique(as.character(obj$sample_id))
sample_ids <- sample_ids[order(sample_ids)]

manifest <- manifest[manifest$sample_id %in% sample_ids, , drop = FALSE]

missing_samples <- setdiff(sample_ids, manifest$sample_id)
if (length(missing_samples) > 0) {
  stop(
    sprintf(
      "The sample manifest is missing sample(s) present in the object: %s",
      paste(missing_samples, collapse = ", ")
    ),
    call. = FALSE
  )
}

sample_objects <- list()
summary_list <- list()

for (sid in sample_ids) {
  log_message("Running SoupX for sample: ", sid)

  sample_obj <- subset(obj, subset = sample_id == sid)
  sample_cells_prefixed <- colnames(sample_obj)
  raw_barcodes <- extract_raw_barcodes(sample_cells_prefixed, sid)

  man_row <- manifest[manifest$sample_id == sid, , drop = FALSE]
  if (nrow(man_row) != 1) {
    stop(sprintf("Expected exactly one manifest row for sample %s", sid), call. = FALSE)
  }

  filtered_counts <- read_10x_flexible(man_row$filtered_path[[1]])
  raw_counts <- read_10x_flexible(man_row$raw_path[[1]])

  keep_barcodes <- intersect(colnames(filtered_counts), raw_barcodes)
  if (length(keep_barcodes) == 0) {
    stop(sprintf("No overlapping filtered barcodes found for sample %s", sid), call. = FALSE)
  }

  sample_clusters <- setNames(as.character(sample_obj$pre_soupx_cluster), raw_barcodes)
  sample_clusters <- sample_clusters[keep_barcodes]
  filtered_counts <- filtered_counts[, keep_barcodes, drop = FALSE]

  if (!identical(names(sample_clusters), colnames(filtered_counts))) {
    stop(sprintf("Cluster labels do not match filtered matrix barcodes for sample %s", sid), call. = FALSE)
  }

  log_message(
    "Sample ", sid,
    ": filtered cells used by SoupX = ", ncol(filtered_counts),
    ", raw droplets = ", ncol(raw_counts)
  )

  sc <- SoupChannel(tod = raw_counts, toc = filtered_counts)
  sc <- setClusters(sc, sample_clusters)

  sc <- tryCatch(
    {
      autoEstCont(sc, doPlot = FALSE)
    },
    error = function(e) {
      log_message("autoEstCont failed for ", sid, "; using fixed contamination fraction 0.10")
      setContaminationFraction(sc, 0.10)
    }
  )

  corrected_counts <- adjustCounts(sc, roundToInt = TRUE)
  if (!inherits(corrected_counts, "dgCMatrix")) {
    corrected_counts <- as(corrected_counts, "dgCMatrix")
  }

  meta_df <- sample_obj@meta.data
  rownames(meta_df) <- raw_barcodes
  meta_df <- meta_df[colnames(corrected_counts), , drop = FALSE]

  if (!identical(rownames(meta_df), colnames(corrected_counts))) {
    stop(sprintf("Metadata rownames do not match corrected count barcodes for sample %s", sid), call. = FALSE)
  }

  corrected_obj <- CreateSeuratObject(
    counts = corrected_counts,
    project = sid,
    meta.data = meta_df
  )

  corrected_obj$sample_id <- sid

  sample_objects[[sid]] <- corrected_obj

  rho_value <- NA_real_
  if (!is.null(sc$fit) && !is.null(sc$fit$rhoEst)) {
    rho_value <- as.numeric(sc$fit$rhoEst)
  }

  summary_list[[sid]] <- data.frame(
    sample_id = sid,
    n_cells = ncol(corrected_obj),
    rho_est = rho_value,
    stringsAsFactors = FALSE
  )

  log_message("Completed SoupX for sample: ", sid)
}

if (length(sample_objects) == 0) {
  stop("No SoupX-corrected sample objects were created.", call. = FALSE)
}

log_message("Merging SoupX-corrected samples")
if (length(sample_objects) == 1) {
  merged_obj <- sample_objects[[1]]
} else {
  merged_obj <- merge(
    x = sample_objects[[1]],
    y = sample_objects[2:length(sample_objects)],
    add.cell.ids = names(sample_objects),
    project = "RetAngio_SoupX"
  )
}

if (!is.null(obj@misc)) {
  merged_obj@misc <- obj@misc
}
merged_obj@misc$soupx <- list(
  input_object = input_rds,
  output_object = output_rds,
  contamination_summary = output_csv,
  cluster_column = "pre_soupx_cluster",
  run_datetime = as.character(Sys.time())
)

summary_df <- do.call(rbind, summary_list)

saveRDS(merged_obj, output_rds)
write.csv(summary_df, output_csv, row.names = FALSE)

log_message("Saved SoupX-corrected object: ", output_rds)
log_message("Saved contamination summary: ", output_csv)
log_message("Stage 05b completed successfully")
