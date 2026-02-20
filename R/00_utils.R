suppressPackageStartupMessages({
  library(yaml)
  library(fs)
  library(glue)
  library(digest)
  library(Seurat)
  library(ggplot2)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) {
    return(list(config = "config/config.yml"))
  }
  key <- which(args == "--config")
  if (length(key) == 1 && key < length(args)) {
    return(list(config = args[key + 1]))
  }
  stop("Usage: Rscript <script> --config config/config.yml", call. = FALSE)
}

load_config <- function(config_path) {
  if (!file.exists(config_path)) {
    stop(glue("Config file not found: {config_path}"), call. = FALSE)
  }
  cfg <- yaml::read_yaml(config_path)
  cfg
}

init_dirs <- function(cfg) {
  out <- cfg$output
  dirs <- file.path(out$root, unlist(out[c("snapshots_dir", "reports_dir", "plots_dir", "tables_dir", "logs_dir")]))
  fs::dir_create(out$root)
  purrr::walk(dirs, fs::dir_create)
  invisible(dirs)
}

log_message <- function(msg, cfg, script_name = NULL) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  prefix <- if (!is.null(script_name)) glue("[{script_name}]") else ""
  line <- glue("{ts} {prefix} {msg}")
  cat(line, "\n")
  log_file <- file.path(cfg$output$root, cfg$output$logs_dir, "pipeline.log")
  cat(line, "\n", file = log_file, append = TRUE)
}

time_it <- function(expr, label, cfg, script_name = NULL) {
  start <- Sys.time()
  on.exit({
    elapsed <- round(as.numeric(difftime(Sys.time(), start, units = "secs")), 2)
    log_message(glue("{label} completed in {elapsed}s"), cfg, script_name)
  })
  force(expr)
}

safe_save_rds <- function(object, path) {
  fs::dir_create(path_dir(path))
  saveRDS(object, path)
  checksum <- digest::digest(file = path, algo = "sha256")
  attr(object, "sha256") <- checksum
  checksum
}

write_session_info <- function(cfg, script_name) {
  out_file <- file.path(cfg$output$root, cfg$output$reports_dir, glue("sessionInfo_{script_name}.txt"))
  sink(out_file)
  print(sessionInfo())
  sink()
}

object_size_report <- function(seu) {
  tibble::tibble(
    object = class(seu)[1],
    cells = ncol(seu),
    features = nrow(seu),
    size_mb = round(as.numeric(object.size(seu)) / 1024^2, 2),
    assays = paste(Assays(seu), collapse = ";"),
    reductions = paste(names(seu@reductions), collapse = ";")
  )
}

validate_seurat <- function(seu, cfg) {
  stopifnot(inherits(seu, "Seurat"))
  required_assays <- cfg$input_assays$required
  missing_assays <- setdiff(required_assays, Assays(seu))
  if (length(missing_assays) > 0) {
    stop(glue("Missing required assays: {paste(missing_assays, collapse = ', ')}"), call. = FALSE)
  }
  cluster_col <- cfg$contaminant_filter$cluster_column
  if (!(cluster_col %in% colnames(seu@meta.data))) {
    stop(glue("Missing cluster column: {cluster_col}"), call. = FALSE)
  }
  invisible(TRUE)
}

compactify_seurat <- function(seu, cfg) {
  # Seurat v5 object compaction by dropping heavy matrices not required downstream.
  # Use carefully; retain RNA data layer for DEG-related workflows.
  assays_keep <- cfg$base_object$assays_to_keep
  for (assay_name in setdiff(Assays(seu), assays_keep)) {
    seu[[assay_name]] <- NULL
  }

  if (isTRUE(cfg$base_object$drop_scale_data)) {
    for (assay_name in Assays(seu)) {
      assay_obj <- seu[[assay_name]]
      if (inherits(assay_obj, "Assay")) {
        assay_obj@scale.data <- matrix(0, nrow = 0, ncol = 0)
        seu[[assay_name]] <- assay_obj
      }
    }
  }
  gc()
  seu
}

save_plot <- function(plot_obj, filename, cfg) {
  out_path <- file.path(cfg$output$root, cfg$output$plots_dir, filename)
  ggplot2::ggsave(
    filename = out_path,
    plot = plot_obj,
    width = cfg$plots$width,
    height = cfg$plots$height,
    dpi = cfg$plots$dpi
  )
  out_path
}

assert_file_nonempty <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop(glue("Expected non-empty file missing: {path}"), call. = FALSE)
  }
}
