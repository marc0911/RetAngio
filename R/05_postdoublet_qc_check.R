source("R/00_utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tibble)
})

extract_qc <- function(seu, sample_col, stage_label) {
  DefaultAssay(seu) <- "RNA"

  if (!"percent.mt" %in% colnames(seu@meta.data)) {
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
  }
  if (!"percent.ribo" %in% colnames(seu@meta.data)) {
    seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^Rpl|^Rps")
  }

  seu@meta.data %>%
    rownames_to_column("cell_id") %>%
    transmute(
      cell_id = cell_id,
      sample_id = .data[[sample_col]],
      stage = stage_label,
      nCount_RNA = nCount_RNA,
      nFeature_RNA = nFeature_RNA,
      percent.mt = percent.mt,
      percent.ribo = percent.ribo
    )
}

main <- function() {
  log_message("Starting stage 05_postdoublet_qc_check")

  cfg <- read_config()
  outdir <- cfg$project$output_dir
  sample_col <- cfg$analysis$sample_column

  ensure_dir(file.path(outdir, "doublets"))
  ensure_dir(file.path("plots", "doublets"))

  pre_rds  <- file.path(outdir, "objects", "02_qc_filtered.rds")
  post_rds <- file.path(outdir, "objects", "04_qc_doublet_filtered.rds")

  if (!file.exists(pre_rds))  stop("Missing input object: ", pre_rds)
  if (!file.exists(post_rds)) stop("Missing input object: ", post_rds)

  seu_pre  <- readRDS(pre_rds)
  seu_post <- readRDS(post_rds)

  log_message("Loaded pre-doublet object: ", pre_rds)
  log_message("Loaded post-doublet object: ", post_rds)

  df_pre  <- extract_qc(seu_pre, sample_col, "pre_doublet_filter")
  df_post <- extract_qc(seu_post, sample_col, "post_doublet_filter")

  df <- bind_rows(df_pre, df_post)

  write_csv(df, file.path(outdir, "doublets", "05_postdoublet_qc_comparison_per_cell.csv.gz"))

  summary_tbl <- df %>%
    group_by(sample_id, stage) %>%
    summarise(
      n_cells = n(),
      median_nCount_RNA = median(nCount_RNA),
      median_nFeature_RNA = median(nFeature_RNA),
      median_percent_mt = median(percent.mt),
      p95_percent_mt = quantile(percent.mt, 0.95),
      .groups = "drop"
    )

  write_csv(summary_tbl, file.path(outdir, "doublets", "05_postdoublet_qc_summary.csv"))

  p_feat <- ggplot(df, aes(x = stage, y = nFeature_RNA)) +
    geom_violin(scale = "width", trim = TRUE) +
    facet_wrap(~ sample_id, scales = "free_y") +
    theme_bw()

  ggsave(
    filename = file.path("plots", "doublets", "05_postdoublet_nFeature_violin.png"),
    plot = p_feat,
    width = 12,
    height = 8,
    dpi = 300
  )

  p_count <- ggplot(df, aes(x = stage, y = nCount_RNA)) +
    geom_violin(scale = "width", trim = TRUE) +
    facet_wrap(~ sample_id, scales = "free_y") +
    theme_bw()

  ggsave(
    filename = file.path("plots", "doublets", "05_postdoublet_nCount_violin.png"),
    plot = p_count,
    width = 12,
    height = 8,
    dpi = 300
  )

  p_mt <- ggplot(df, aes(x = stage, y = percent.mt)) +
    geom_violin(scale = "width", trim = TRUE) +
    facet_wrap(~ sample_id, scales = "free_y") +
    theme_bw()

  ggsave(
    filename = file.path("plots", "doublets", "05_postdoublet_percent_mt_violin.png"),
    plot = p_mt,
    width = 12,
    height = 8,
    dpi = 300
  )

  p_cells <- summary_tbl %>%
    ggplot(aes(x = sample_id, y = n_cells, fill = stage)) +
    geom_col(position = "dodge") +
    theme_bw()

  ggsave(
    filename = file.path("plots", "doublets", "05_postdoublet_cell_counts.png"),
    plot = p_cells,
    width = 9,
    height = 5,
    dpi = 300
  )

  log_message("Stage 05_postdoublet_qc_check complete")
}

main()
