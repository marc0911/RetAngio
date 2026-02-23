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

get_io_config <- function(cfg) {
  io <- cfg$io
  if (is.null(io)) {
    io <- list(backend = "rds", rds_compress = "gzip", fst_compress = 50)
  }

  io$backend <- io$backend %||% "rds"
  io$rds_compress <- io$rds_compress %||% "gzip"
  io$fst_compress <- io$fst_compress %||% 50

  allowed_backends <- c("rds", "qs2", "fst_df")
  if (!io$backend %in% allowed_backends) {
    fail_fast("Invalid io.backend={io$backend}; allowed: {paste(allowed_backends, collapse=', ')}")
  }

  valid_rds <- is.logical(io$rds_compress) || io$rds_compress %in% c("gzip", "bzip2", "xz")
  if (!valid_rds) {
    fail_fast("io.rds_compress must be TRUE/FALSE or one of: gzip, bzip2, xz")
  }

  io
}

serialization_extension <- function(backend) {
  dplyr::case_when(
    backend == "rds" ~ "rds",
    backend == "qs2" ~ "qs2",
    backend == "fst_df" ~ "fst",
    TRUE ~ NA_character_
  )
}

serialization_path <- function(stem, cfg) {
  io <- get_io_config(cfg)
  ext <- serialization_extension(io$backend)
  stem <- sub("\\.(rds|qs2|fst)$", "", stem)
  paste0(stem, ".", ext)
}

save_object <- function(object, stem, cfg, object_type = c("auto", "seurat", "data.frame", "generic")) {
  io <- get_io_config(cfg)
  object_type <- match.arg(object_type)

  if (object_type == "auto") {
    object_type <- if (inherits(object, "Seurat")) "seurat" else if (is.data.frame(object)) "data.frame" else "generic"
  }

  backend <- io$backend
  if (backend == "fst_df" && object_type != "data.frame") {
    warning("io.backend=fst_df supports data.frame only; falling back to RDS for non-data.frame object", call. = FALSE)
    backend <- "rds"
  }
  if (backend == "qs2" && !requireNamespace("qs2", quietly = TRUE)) {
    warning("io.backend=qs2 requested but package 'qs2' is not installed; falling back to RDS", call. = FALSE)
    backend <- "rds"
  }
  if (backend == "fst_df" && !requireNamespace("fst", quietly = TRUE)) {
    warning("io.backend=fst_df requested but package 'fst' is not installed; falling back to RDS", call. = FALSE)
    backend <- "rds"
  }

  path <- paste0(sub("\\.(rds|qs2|fst)$", "", stem), ".", serialization_extension(backend))

  if (backend == "rds") {
    # Base R serialization documentation: ?saveRDS / ?readRDS
    saveRDS(object, file = path, compress = io$rds_compress)
  } else if (backend == "qs2") {
    # qs2 package documentation (CRAN): ?qs2::qs_save / ?qs2::qs_read
    qs2::qs_save(object, file = path)
  } else if (backend == "fst_df") {
    # fst documentation: write_fst/read_fst are for tabular data frames.
    fst::write_fst(object, path = path, compress = io$fst_compress)
  }

  log_info("Saved object with backend={backend} at {path}")
  path
}

load_object <- function(stem_or_path, cfg = NULL) {
  if (!is.null(cfg)) {
    io <- get_io_config(cfg)
    backend <- io$backend
    if (backend == "qs2" && !requireNamespace("qs2", quietly = TRUE)) {
      warning("io.backend=qs2 requested but package 'qs2' is not installed; falling back to RDS", call. = FALSE)
      backend <- "rds"
    }
    if (backend == "fst_df" && !requireNamespace("fst", quietly = TRUE)) {
      warning("io.backend=fst_df requested but package 'fst' is not installed; falling back to RDS", call. = FALSE)
      backend <- "rds"
    }
    path <- paste0(sub("\\.(rds|qs2|fst)$", "", stem_or_path), ".", serialization_extension(backend))
  } else {
    path <- stem_or_path
    ext <- tools::file_ext(path)
    backend <- dplyr::case_when(
      ext == "rds" ~ "rds",
      ext == "qs2" ~ "qs2",
      ext == "fst" ~ "fst_df",
      TRUE ~ "unknown"
    )
    if (backend == "unknown") {
      fail_fast("Cannot infer backend from extension: {path}")
    }
  }

  if (!file.exists(path)) {
    fail_fast("Serialized file not found: {path}")
  }

  obj <- if (backend == "rds") {
    readRDS(path)
  } else if (backend == "qs2") {
    if (!requireNamespace("qs2", quietly = TRUE)) {
      fail_fast("Cannot load qs2 file; package 'qs2' is not installed: {path}")
    }
    qs2::qs_read(path)
  } else if (backend == "fst_df") {
    if (!requireNamespace("fst", quietly = TRUE)) {
      fail_fast("Cannot load fst file; package 'fst' is not installed: {path}")
    }
    fst::read_fst(path, as.data.table = FALSE)
  }

  log_info("Loaded object with backend={backend} from {path}")
  obj
}

run_io_sanity_checks <- function(cfg, outdir) {
  ensure_dir(outdir)

  small_list <- list(a = 1L, b = "retangio", c = c(0.1, 0.2))
  small_df <- tibble::tibble(x = 1:3, y = c("a", "b", "c"))

  backends <- c("rds")
  if (requireNamespace("qs2", quietly = TRUE)) backends <- c(backends, "qs2")
  if (requireNamespace("fst", quietly = TRUE)) backends <- c(backends, "fst_df")

  rows <- list()
  for (b in backends) {
    cfg_b <- cfg
    cfg_b$io$backend <- b

    list_path <- save_object(small_list, file.path(outdir, paste0("00_io_sanity_list_", b)), cfg_b, object_type = "generic")
    list_back <- load_object(list_path)
    ok_list <- isTRUE(all.equal(small_list, list_back))

    df_path <- save_object(small_df, file.path(outdir, paste0("00_io_sanity_df_", b)), cfg_b, object_type = "data.frame")
    df_back <- load_object(df_path)
    ok_df <- isTRUE(all.equal(as.data.frame(small_df), as.data.frame(df_back)))

    rows[[length(rows) + 1]] <- tibble::tibble(
      backend = b,
      list_path = list_path,
      list_equal = ok_list,
      df_path = df_path,
      df_equal = ok_df
    )
    log_info("IO sanity backend={b}; list_equal={ok_list}; df_equal={ok_df}")
  }

  summary <- dplyr::bind_rows(rows)
  readr::write_csv(summary, file.path(outdir, "00_io_sanity_checks.csv"))

  if (!all(summary$list_equal, summary$df_equal)) {
    fail_fast("IO sanity checks failed; see outputs/00_io_sanity_checks.csv")
  }

  invisible(summary)
}


`%||%` <- function(x, y) if (is.null(x)) y else x

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
