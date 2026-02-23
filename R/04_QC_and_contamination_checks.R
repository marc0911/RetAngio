source("R/00_utils.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
})

safe_feature <- function(seu, genes, assay = "RNA") {
  present <- intersect(genes, rownames(seu[[assay]]))
  if (length(present) == 0) return(rep(0, ncol(seu)))
  Matrix::colSums(get_layer_data(seu, assay, "data", fallback_slot = "data")[present, , drop = FALSE])
}

main <- function() {
  cfg <- read_config()
  outdir <- cfg$project$output_dir
  ensure_dir(outdir)

  seu <- load_object(file.path(outdir, "03_BASE_object"), cfg)
  check_required_layers(seu, "RNA")

  seu$score_hb <- safe_feature(seu, cfg$qc$hemoglobin_genes)
  seu$score_photo <- safe_feature(seu, cfg$qc$photoreceptor_markers)
  seu$score_rpe <- safe_feature(seu, cfg$qc$rpe_markers)

  qc_tbl <- seu@meta.data |>
    tibble::rownames_to_column("cell") |>
    group_by(final_cluster) |>
    summarise(
      n_cells = n(),
      pct_mt_mean = mean(percent.mt, na.rm = TRUE),
      pct_ribo_mean = mean(percent.ribo, na.rm = TRUE),
      hb_score_mean = mean(score_hb, na.rm = TRUE),
      photo_score_mean = mean(score_photo, na.rm = TRUE),
      rpe_score_mean = mean(score_rpe, na.rm = TRUE),
      .groups = "drop"
    )
  write_csv(qc_tbl, file.path(outdir, "04_qc_cluster_summary.csv"))

  p1 <- VlnPlot(seu, features = c("percent.mt", "percent.ribo"), group.by = "final_cluster", pt.size = 0)
  p2 <- VlnPlot(seu, features = c("score_hb", "score_photo", "score_rpe"), group.by = "final_cluster", pt.size = 0)
  ggsave(file.path(outdir, "04_qc_violins.png"), p1 / p2, width = 12, height = 10, dpi = 300)

  # Optional SoupX branch (Young & Behjati 2020, PMID:33367645)
  if (isTRUE(cfg$soupx$enabled)) {
    warning("SoupX requested but this pipeline only logs requirement: provide raw matrices in cfg$soupx$raw_mtx_dir", call. = FALSE)
  }

  log_info("04_QC_and_contamination_checks complete")
}

if (sys.nframe() == 0) {
  main()
}
