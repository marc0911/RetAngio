source("R/00_utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tibble)
  library(tidyr)
})

main <- function() {
  log_message("Starting stage 03_detect_doublets")

  cfg <- read_config()
  outdir <- cfg$project$output_dir
  sample_col <- cfg$analysis$sample_column

  ensure_dir(outdir)
  ensure_dir(file.path(outdir, "objects"))
  ensure_dir(file.path(outdir, "doublets"))
  ensure_dir("plots")
  ensure_dir(file.path("plots", "doublets"))

  in_rds <- file.path(outdir, "objects", "02_qc_filtered.rds")
  if (!file.exists(in_rds)) {
    stop("Missing input object: ", in_rds)
  }

  seu <- readRDS(in_rds)
  log_message("Loaded object: ", in_rds)

  if (!sample_col %in% colnames(seu@meta.data)) {
    stop("Sample column not found in metadata: ", sample_col)
  }

  DefaultAssay(seu) <- "RNA"

  # Seurat v5 safeguard: join RNA layers before conversion
  if ("RNA" %in% names(seu@assays)) {
    log_message("Joining RNA layers before SingleCellExperiment conversion")
    seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  }

  sce <- as.SingleCellExperiment(seu, assay = "RNA")
  log_message("Converted Seurat -> SingleCellExperiment")

  set.seed(cfg$project$seed)

  sce <- scDblFinder(
    sce,
    samples = colData(sce)[[sample_col]],
    verbose = TRUE
  )

  colData_df <- as.data.frame(colData(sce))

  needed_cols <- c("scDblFinder.score", "scDblFinder.class")
  missing_cols <- setdiff(needed_cols, colnames(colData_df))
  if (length(missing_cols) > 0) {
    stop("Missing expected scDblFinder columns: ", paste(missing_cols, collapse = ", "))
  }

  seu$scDblFinder.score <- colData_df$scDblFinder.score
  seu$scDblFinder.class <- colData_df$scDblFinder.class

  saveRDS(seu, file.path(outdir, "objects", "03_qc_filtered_doublets_annotated.rds"))
  log_message("Saved object: ", file.path(outdir, "objects", "03_qc_filtered_doublets_annotated.rds"))

  meta <- seu@meta.data %>%
    rownames_to_column("cell_id") %>%
    transmute(
      cell_id = cell_id,
      sample_id = .data[[sample_col]],
      scDblFinder.score = scDblFinder.score,
      scDblFinder.class = scDblFinder.class
    )

  write_csv(meta, file.path(outdir, "doublets", "03_scDblFinder_per_cell.csv.gz"))

  summary_tbl <- meta %>%
    group_by(sample_id, scDblFinder.class) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = scDblFinder.class, values_from = n, values_fill = 0)

  if (!"doublet" %in% colnames(summary_tbl)) summary_tbl$doublet <- 0
  if (!"singlet" %in% colnames(summary_tbl)) summary_tbl$singlet <- 0

  summary_tbl <- summary_tbl %>%
    mutate(
      total = doublet + singlet,
      pct_doublet = 100 * doublet / total
    ) %>%
    arrange(sample_id)

  write_csv(summary_tbl, file.path(outdir, "doublets", "03_scDblFinder_summary_by_sample.csv"))
  log_message("Saved doublet summary by sample")

  p_bar <- ggplot(summary_tbl, aes(x = sample_id, y = pct_doublet)) +
    geom_col() +
    theme_bw() +
    ylab("Percent doublets") +
    xlab("Sample")

  ggsave(
    filename = file.path("plots", "doublets", "03_scDblFinder_pct_doublets_by_sample.png"),
    plot = p_bar,
    width = 8,
    height = 5,
    dpi = 300
  )

  p_hist <- ggplot(meta, aes(x = scDblFinder.score, fill = scDblFinder.class)) +
    geom_histogram(bins = 100, position = "identity", alpha = 0.6) +
    facet_wrap(~ sample_id, scales = "free_y") +
    theme_bw()

  ggsave(
    filename = file.path("plots", "doublets", "03_scDblFinder_score_histograms.png"),
    plot = p_hist,
    width = 12,
    height = 8,
    dpi = 300
  )

  log_message("Stage 03_detect_doublets complete")
}

main()
