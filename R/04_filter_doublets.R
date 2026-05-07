source("R/00_utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tibble)
})

main <- function() {
  log_message("Starting stage 04_filter_doublets")

  cfg <- read_config()
  outdir <- cfg$project$output_dir
  sample_col <- cfg$analysis$sample_column

  ensure_dir(outdir)
  ensure_dir(file.path(outdir, "objects"))
  ensure_dir(file.path(outdir, "doublets"))
  ensure_dir("plots")
  ensure_dir(file.path("plots", "doublets"))

  in_rds <- file.path(outdir, "objects", "03_qc_filtered_doublets_annotated.rds")
  if (!file.exists(in_rds)) {
    stop("Missing input object: ", in_rds)
  }

  seu <- readRDS(in_rds)
  log_message("Loaded object: ", in_rds)

  if (!sample_col %in% colnames(seu@meta.data)) {
    stop("Sample column not found in metadata: ", sample_col)
  }
  if (!"scDblFinder.class" %in% colnames(seu@meta.data)) {
    stop("Metadata column missing: scDblFinder.class")
  }

  meta_before <- seu@meta.data %>%
    rownames_to_column("cell_id") %>%
    transmute(
      cell_id = cell_id,
      sample_id = .data[[sample_col]],
      scDblFinder.class = scDblFinder.class
    )

  seu_filt <- subset(seu, subset = scDblFinder.class == "singlet")

  meta_after <- seu_filt@meta.data %>%
    rownames_to_column("cell_id") %>%
    transmute(
      cell_id = cell_id,
      sample_id = .data[[sample_col]]
    )

  summary_before <- meta_before %>%
    count(sample_id, name = "n_before")

  summary_after <- meta_after %>%
    count(sample_id, name = "n_after")

  doublet_tbl <- meta_before %>%
    count(sample_id, scDblFinder.class, name = "n") %>%
    tidyr::pivot_wider(names_from = scDblFinder.class, values_from = n, values_fill = 0)

  if (!"doublet" %in% colnames(doublet_tbl)) doublet_tbl$doublet <- 0
  if (!"singlet" %in% colnames(doublet_tbl)) doublet_tbl$singlet <- 0

  summary_tbl <- summary_before %>%
    left_join(summary_after, by = "sample_id") %>%
    left_join(doublet_tbl, by = "sample_id") %>%
    mutate(
      n_after = ifelse(is.na(n_after), 0L, n_after),
      n_removed = n_before - n_after,
      pct_removed = ifelse(n_before > 0, 100 * n_removed / n_before, NA_real_)
    ) %>%
    arrange(sample_id)

  saveRDS(seu_filt, file.path(outdir, "objects", "04_qc_doublet_filtered.rds"))
  log_message("Saved object: ", file.path(outdir, "objects", "04_qc_doublet_filtered.rds"))

  write_csv(summary_tbl, file.path(outdir, "doublets", "04_doublet_filtering_summary.csv"))

  p_bar <- ggplot(summary_tbl, aes(x = sample_id, y = pct_removed)) +
    geom_col() +
    theme_bw() +
    ylab("Percent removed as doublets") +
    xlab("Sample")

  ggsave(
    filename = file.path("plots", "doublets", "04_pct_removed_doublets_by_sample.png"),
    plot = p_bar,
    width = 8,
    height = 5,
    dpi = 300
  )

  log_message("Cells before doublet filtering: ", ncol(seu))
  log_message("Cells after doublet filtering: ", ncol(seu_filt))
  log_message("Stage 04_filter_doublets complete")
}

main()
