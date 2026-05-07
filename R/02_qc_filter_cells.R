source("R/00_utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(patchwork)
  library(yaml)
})

main <- function() {
  log_message("Starting stage 02_qc_filter_cells")

  cfg <- read_config()
  outdir <- cfg$project$output_dir

  ensure_dir(outdir)
  ensure_dir(file.path(outdir, "objects"))
  ensure_dir(file.path(outdir, "qc"))
  ensure_dir("plots")
  ensure_dir(file.path("plots", "qc"))

  rds_path <- file.path(outdir, "objects", "00_raw_merged.rds")
  decision_path <- "decisions/01_qc_thresholds.yml"

  if (!file.exists(rds_path)) {
    stop("Missing input object: ", rds_path)
  }
  if (!file.exists(decision_path)) {
    stop("Missing decision file: ", decision_path)
  }

  seu <- readRDS(rds_path)
  log_message("Loaded object: ", rds_path)

  dec <- yaml::read_yaml(decision_path)
  sample_col <- cfg$analysis$sample_column

  if (!sample_col %in% colnames(seu@meta.data)) {
    stop("Sample column not found in metadata: ", sample_col)
  }

  DefaultAssay(seu) <- "RNA"

  # Recompute if missing
  if (!"percent.mt" %in% colnames(seu@meta.data)) {
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
  }
  if (!"percent.ribo" %in% colnames(seu@meta.data)) {
    seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^Rpl|^Rps")
  }

  min_features <- dec$decision$min_features_rna
  max_percent_mt <- dec$decision$max_percent_mt
  min_counts <- dec$decision$min_counts_rna
  max_counts <- dec$decision$max_counts_rna
  max_features <- dec$decision$max_features_rna

  meta_before <- seu@meta.data %>%
    rownames_to_column("cell_id") %>%
    transmute(
      cell_id = cell_id,
      sample_id = .data[[sample_col]],
      nCount_RNA = nCount_RNA,
      nFeature_RNA = nFeature_RNA,
      percent.mt = percent.mt,
      percent.ribo = percent.ribo
    )

  keep <- rep(TRUE, ncol(seu))

  if (!is.null(min_features)) {
    keep <- keep & (seu$nFeature_RNA >= min_features)
  }
  if (!is.null(max_features)) {
    keep <- keep & (seu$nFeature_RNA <= max_features)
  }
  if (!is.null(min_counts)) {
    keep <- keep & (seu$nCount_RNA >= min_counts)
  }
  if (!is.null(max_counts)) {
    keep <- keep & (seu$nCount_RNA <= max_counts)
  }
  if (!is.null(max_percent_mt)) {
    keep <- keep & (seu$percent.mt <= max_percent_mt)
  }

  kept_cells <- colnames(seu)[keep]
  seu_filt <- subset(seu, cells = kept_cells)

  log_message("Cells before filtering: ", ncol(seu))
  log_message("Cells after filtering: ", ncol(seu_filt))

  meta_after <- seu_filt@meta.data %>%
    rownames_to_column("cell_id") %>%
    transmute(
      cell_id = cell_id,
      sample_id = .data[[sample_col]],
      nCount_RNA = nCount_RNA,
      nFeature_RNA = nFeature_RNA,
      percent.mt = percent.mt,
      percent.ribo = percent.ribo
    )

  summary_before <- meta_before %>%
    count(sample_id, name = "n_before")

  summary_after <- meta_after %>%
    count(sample_id, name = "n_after")

  summary_tbl <- summary_before %>%
    full_join(summary_after, by = "sample_id") %>%
    mutate(
      n_before = ifelse(is.na(n_before), 0L, n_before),
      n_after = ifelse(is.na(n_after), 0L, n_after),
      n_removed = n_before - n_after,
      pct_removed = ifelse(n_before > 0, 100 * n_removed / n_before, NA_real_)
    ) %>%
    arrange(sample_id)

  write_csv(summary_tbl, file.path(outdir, "qc", "02_qc_filtering_summary.csv"))
  write_csv(meta_after, file.path(outdir, "qc", "02_qc_filtered_cell_metrics.csv.gz"))

  saveRDS(seu_filt, file.path(outdir, "objects", "02_qc_filtered.rds"))
  log_message("Saved filtered object: ", file.path(outdir, "objects", "02_qc_filtered.rds"))

  # Post-filter QC violin plots
  p_vln1 <- VlnPlot(
    seu_filt,
    features = c("nFeature_RNA", "nCount_RNA"),
    group.by = sample_col,
    pt.size = 0
  ) + patchwork::plot_layout(ncol = 2)

  ggsave(
    filename = file.path("plots", "qc", "02_postfilter_violin_features.png"),
    plot = p_vln1,
    width = 12,
    height = 6,
    dpi = 300
  )

  p_vln2 <- VlnPlot(
    seu_filt,
    features = c("percent.mt", "percent.ribo"),
    group.by = sample_col,
    pt.size = 0
  ) + patchwork::plot_layout(ncol = 2)

  ggsave(
    filename = file.path("plots", "qc", "02_postfilter_violin_percentages.png"),
    plot = p_vln2,
    width = 12,
    height = 6,
    dpi = 300
  )

  # Filtering barplot
  p_bar <- ggplot(summary_tbl, aes(x = sample_id, y = n_removed)) +
    geom_col() +
    theme_bw() +
    ylab("Removed cells") +
    xlab("Sample")

  ggsave(
    filename = file.path("plots", "qc", "02_qc_removed_cells_by_sample.png"),
    plot = p_bar,
    width = 8,
    height = 5,
    dpi = 300
  )

  log_message("Stage 02_qc_filter_cells complete")
}

main()
