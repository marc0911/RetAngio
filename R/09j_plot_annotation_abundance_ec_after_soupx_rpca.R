#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

source("R/plotting_helpers.R")

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_rds <- "results/objects/09h_ec_after_soupx_rpca_annotated.rds"
  outdir_plots <- "plots/09j_ec_after_soupx_rpca_annotation_abundance"
  outdir_results <- "results/abundance"

  dir.create(outdir_plots, recursive = TRUE, showWarnings = FALSE)
  dir.create(outdir_results, recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 09j_plot_annotation_abundance_ec_after_soupx_rpca")

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  seu <- retangio_add_plotting_metadata(seu)

  if (!"annotation_my" %in% colnames(seu@meta.data)) {
    stop("annotation_my not found in metadata.")
  }

  # NOTE:
  # I interpreted the duplicated 'Arterial x capillary' from the earlier message
  # as a typo and used 'Uncertain' as the final class in the order below.
  label_order <- c(
    "Arterial",
    "Arterial x capillary",
    "Capillary 1",
    "Capillary 2",
    "Venous x capillary",
    "Tip cell",
    "Proliferative 1",
    "Proliferative 2",
    "Proliferative 3",
    "Inflammatory EC",
    "EC1",
    "EC2",
    "Uncertain"
  )

  seu$annotation_my <- factor(as.character(seu$annotation_my), levels = label_order)

  palette_values <- retangio_make_named_palette(label_order)

  retangio_write_annotation_abundance_outputs(
    obj = seu,
    annotation_col = "annotation_my",
    label_order = label_order,
    outdir_plots = outdir_plots,
    outdir_results = outdir_results,
    plot_prefix = "09j_rpca",
    table_prefix = "09j_ec_after_soupx_rpca",
    palette_values = palette_values,
    title_prefix = "RPCA — annotated EC subtypes"
  )

  message_ts("Stage 09j_plot_annotation_abundance_ec_after_soupx_rpca complete")
}

main()
