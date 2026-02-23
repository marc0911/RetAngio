suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(glue)
  library(readr)
  library(digest)
  library(Matrix)
  library(Seurat)
})

log_info <- function(...) {
  msg <- glue(...)
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg)
}

fail_fast <- function(...) {
  msg <- glue(...)
  stop(msg, call. = FALSE)
}

read_config <- function(path = "config/config.yml") {
  if (!file.exists(path)) {
    fail_fast("Missing config file: {path}")
  }
  yaml::read_yaml(path)
}

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

set_seed_from_config <- function(cfg) {
  seed <- cfg$project$seed
  if (is.null(seed) || !is.numeric(seed)) {
    fail_fast("config project.seed must be numeric")
  }
  set.seed(seed)
  log_info("Using deterministic seed: {seed}")
}

save_config_snapshot <- function(cfg, outdir) {
  ensure_dir(outdir)
  yaml::write_yaml(cfg, file.path(outdir, "config_used.yml"))
}

save_session_info <- function(outdir) {
  ensure_dir(outdir)
  capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))
}

save_package_versions <- function(outdir) {
  ensure_dir(outdir)
  pkgs <- installed.packages()[, c("Package", "Version")]
  readr::write_csv(as.data.frame(pkgs), file.path(outdir, "package_versions.csv"))
}

get_layer_data <- function(seu, assay, layer, fallback_slot = NULL) {
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    return(SeuratObject::LayerData(seu[[assay]], layer = layer))
  }
  if (is.null(fallback_slot)) {
    fail_fast("SeuratObject<5 requires fallback_slot for assay={assay} layer={layer}")
  }
  warning(glue("Using slot={fallback_slot} fallback for SeuratObject<5"), call. = FALSE)
  SeuratObject::GetAssayData(seu[[assay]], slot = fallback_slot)
}

set_layer_data <- function(seu, assay, layer, value, fallback_slot = NULL) {
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    SeuratObject::LayerData(seu[[assay]], layer = layer) <- value
    return(seu)
  }
  if (is.null(fallback_slot)) {
    fail_fast("SeuratObject<5 requires fallback_slot for assay={assay} layer={layer}")
  }
  warning(glue("Using slot={fallback_slot} fallback for SeuratObject<5"), call. = FALSE)
  SeuratObject::SetAssayData(seu, assay = assay, slot = fallback_slot, new.data = value)
}

check_required_layers <- function(seu, assay = "RNA") {
  layers <- SeuratObject::Layers(seu[[assay]])
  required <- c("counts", "data")
  missing_layers <- setdiff(required, layers)
  if (length(missing_layers) > 0) {
    fail_fast("Assay {assay} missing required layers: {paste(missing_layers, collapse = ', ')}")
  }
}

write_checksum <- function(path, out_file) {
  tibble::tibble(file = basename(path), md5 = digest::digest(file = path, algo = "md5")) |>
    readr::write_csv(out_file)
}

check_contaminants_removed <- function(seu, cluster_col, contaminants) {
  cl <- as.character(seu[[cluster_col]][, 1])
  if (any(cl %in% as.character(contaminants))) {
    fail_fast("Contaminants still present in {cluster_col}: {paste(unique(cl[cl %in% contaminants]), collapse=', ')}")
  }
  TRUE
}

assert_ident_matches <- function(seu, cluster_col) {
  idents <- as.character(Seurat::Idents(seu))
  cols <- as.character(seu[[cluster_col]][, 1])
  if (!identical(idents, cols)) {
    fail_fast("Idents do not match {cluster_col}")
  }
  TRUE
}

assert_reductions_present <- function(seu, reduction_names = c("pca", "integrated.rpca", "umap")) {
  missing <- setdiff(reduction_names, names(seu@reductions))
  if (length(missing) > 0) {
    fail_fast("Missing reductions: {paste(missing, collapse=', ')}")
  }
  TRUE
}
