source("R/00_utils.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(purrr)
  library(tidyr)
})

# Method notes:
# - Ambient RNA correction uses SoupX when raw droplet matrices are supplied.
#   Reference: Young & Behjati 2020, GigaScience, DOI:10.1093/gigascience/giaa151.

args <- parse_args()
cfg <- load_config(args$config)
set.seed(cfg$seed)
init_dirs(cfg)
script_name <- "04_QC_and_contamination_checks"

infile <- file.path(cfg$output$root, cfg$output$snapshots_dir, "03_BASE_object.rds")
if (!file.exists(infile)) stop("Missing 03 output.", call. = FALSE)
seu <- readRDS(infile)

cluster_col <- paste0("SCT_snn_res.", cfg$recluster$final_resolution)
if (!(cluster_col %in% colnames(seu@meta.data))) {
  cluster_col <- tail(grep("^SCT_snn_res\\.", colnames(seu@meta.data), value = TRUE), 1)
}

qc_features <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")
for (feature in qc_features) {
  if (feature %in% colnames(seu@meta.data)) {
    p <- VlnPlot(seu, features = feature, group.by = "orig.ident", pt.size = 0) + ggtitle(paste("QC by sample:", feature))
    save_plot(p, paste0("04_qc_sample_", feature, ".png"), cfg)

    p2 <- VlnPlot(seu, features = feature, group.by = cluster_col, pt.size = 0) + ggtitle(paste("QC by cluster:", feature))
    save_plot(p2, paste0("04_qc_cluster_", feature, ".png"), cfg)
  }
}

contamination_sets <- list(
  hb = cfg$qc$hb_genes,
  photoreceptor = cfg$qc$photoreceptor_markers,
  rpe = cfg$qc$rpe_markers
)

DefaultAssay(seu) <- "RNA"
contam_summary <- purrr::imap_dfr(contamination_sets, function(genes, label) {
  genes_present <- intersect(genes, rownames(seu))
  if (length(genes_present) == 0) {
    return(tibble::tibble(panel = label, gene = NA_character_, pct_expressing = NA_real_, avg_expr = NA_real_))
  }
  expr <- FetchData(seu, vars = genes_present)
  tibble::tibble(
    panel = label,
    gene = genes_present,
    pct_expressing = colMeans(expr > 0),
    avg_expr = colMeans(expr)
  )
})
readr::write_csv(contam_summary, file.path(cfg$output$root, cfg$output$tables_dir, "04_contamination_summary.csv"))

for (panel in names(contamination_sets)) {
  genes <- intersect(contamination_sets[[panel]], rownames(seu))
  if (length(genes) > 0) {
    p <- DotPlot(seu, features = genes, group.by = cluster_col) + RotatedAxis() + ggtitle(paste("Contamination panel:", panel))
    save_plot(p, paste0("04_contam_dotplot_", panel, ".png"), cfg)
  }
}

raw_dirs <- cfg$paths$raw_matrix_dirs
soup_status <- "skipped_no_raw_matrices"
if (length(raw_dirs) > 0) {
  if (!requireNamespace("SoupX", quietly = TRUE)) {
    warning("SoupX not installed; skipping ambient correction.")
    soup_status <- "skipped_missing_package"
  } else {
    # Placeholder implementation: sample-wise loading and correction should be adapted
    # to match provided CellRanger folder structure.
    soup_status <- "configured_but_requires_sample_specific_loading"
  }
}

readr::write_csv(tibble::tibble(step = "soupx", status = soup_status), file.path(cfg$output$root, cfg$output$reports_dir, "04_soupx_status.csv"))

out <- file.path(cfg$output$root, cfg$output$snapshots_dir, "04_qc_checked.rds")
safe_save_rds(seu, out)
write_session_info(cfg, script_name)
log_message("Completed QC + contamination checks", cfg, script_name)
