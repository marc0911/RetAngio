source("R/00_utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(patchwork)
  library(tibble)
})

main <- function() {
  log_message("Starting stage 01_qc_metrics")

  cfg <- read_config()
  outdir <- cfg$project$output_dir

  ensure_dir(outdir)
  ensure_dir(file.path(outdir, "qc"))
  ensure_dir("plots")
  ensure_dir(file.path("plots", "qc"))

  rds_path <- file.path(outdir, "objects", "00_raw_merged.rds")
  if (!file.exists(rds_path)) {
    stop("Missing input object: ", rds_path)
  }

  seu <- readRDS(rds_path)
  log_message("Loaded object: ", rds_path)

  sample_col <- cfg$analysis$sample_column
  if (!sample_col %in% colnames(seu@meta.data)) {
    stop("Sample column not found in metadata: ", sample_col)
  }

  DefaultAssay(seu) <- "RNA"

  # QC metrics
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
  seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^Rpl|^Rps")

  meta <- seu@meta.data %>%
    rownames_to_column("cell_id") %>%
    transmute(
      cell_id = cell_id,
      sample_id = .data[[sample_col]],
      nCount_RNA = nCount_RNA,
      nFeature_RNA = nFeature_RNA,
      percent.mt = percent.mt,
      percent.ribo = percent.ribo
    )

  write_csv(meta, file.path(outdir, "qc", "01_qc_cell_metrics.csv.gz"))
  log_message("Saved per-cell QC metrics")

  summary_tbl <- meta %>%
    group_by(sample_id) %>%
    summarise(
      n_cells = n(),
      median_nCount_RNA = median(nCount_RNA),
      median_nFeature_RNA = median(nFeature_RNA),
      median_percent_mt = median(percent.mt),
      median_percent_ribo = median(percent.ribo),
      p05_nFeature_RNA = quantile(nFeature_RNA, 0.05),
      p95_nFeature_RNA = quantile(nFeature_RNA, 0.95),
      p95_percent_mt = quantile(percent.mt, 0.95),
      .groups = "drop"
    )

  write_csv(summary_tbl, file.path(outdir, "qc", "01_qc_summary_by_sample.csv"))
  log_message("Saved QC summary by sample")

  # Violin plots
  p_vln1 <- VlnPlot(
    seu,
    features = c("nFeature_RNA", "nCount_RNA"),
    group.by = sample_col,
    pt.size = 0
  ) + patchwork::plot_layout(ncol = 2)

  ggsave(
    filename = file.path("plots", "qc", "01_qc_violin_features.png"),
    plot = p_vln1,
    width = 12,
    height = 6,
    dpi = 300
  )

  p_vln2 <- VlnPlot(
    seu,
    features = c("percent.mt", "percent.ribo"),
    group.by = sample_col,
    pt.size = 0
  ) + patchwork::plot_layout(ncol = 2)

  ggsave(
    filename = file.path("plots", "qc", "01_qc_violin_percentages.png"),
    plot = p_vln2,
    width = 12,
    height = 6,
    dpi = 300
  )

  # Scatter plots
  p_sc1 <- FeatureScatter(
    seu,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA",
    group.by = sample_col
  )

  ggsave(
    filename = file.path("plots", "qc", "01_qc_scatter_nCount_vs_nFeature.png"),
    plot = p_sc1,
    width = 7,
    height = 6,
    dpi = 300
  )

  p_sc2 <- FeatureScatter(
    seu,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt",
    group.by = sample_col
  )

  ggsave(
    filename = file.path("plots", "qc", "01_qc_scatter_nCount_vs_percent_mt.png"),
    plot = p_sc2,
    width = 7,
    height = 6,
    dpi = 300
  )

  # Histograms
  p_hist_feat <- ggplot(meta, aes(x = nFeature_RNA)) +
    geom_histogram(bins = 100) +
    facet_wrap(~ sample_id, scales = "free_y") +
    theme_bw()

  ggsave(
    filename = file.path("plots", "qc", "01_qc_hist_nFeature.png"),
    plot = p_hist_feat,
    width = 12,
    height = 8,
    dpi = 300
  )

  p_hist_mt <- ggplot(meta, aes(x = percent.mt)) +
    geom_histogram(bins = 100) +
    facet_wrap(~ sample_id, scales = "free_y") +
    theme_bw()

  ggsave(
    filename = file.path("plots", "qc", "01_qc_hist_percent_mt.png"),
    plot = p_hist_mt,
    width = 12,
    height = 8,
    dpi = 300
  )

  log_message("Stage 01_qc_metrics complete")
}

main()
