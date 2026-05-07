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
  input_rds <- "results/objects/10h_ec_after_soupx_harmony_annotated.rds"
  outdir_plots <- "plots/10j_ec_after_soupx_harmony_annotation_abundance"
  outdir_results <- "results/abundance"

  dir.create(outdir_plots, recursive = TRUE, showWarnings = FALSE)
  dir.create(outdir_results, recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 10j_plot_annotation_abundance_ec_after_soupx_harmony")

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  seu <- retangio_add_plotting_metadata(seu)

  if (!"annotation_my" %in% colnames(seu@meta.data)) {
    stop("annotation_my not found in metadata.")
  }

  label_order <- c(
    "Arterial x Capillary",
    "Capillary",
    "Venous x Capillary",
    "Tip cell",
    "Proliferative 1",
    "Proliferative 2",
    "Inflammatory",
    "EC 1",
    "EC 2",
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
    plot_prefix = "10j_harmony",
    table_prefix = "10j_ec_after_soupx_harmony",
    palette_values = palette_values,
    title_prefix = "Harmony — annotated EC subtypes"
  )

  message_ts("Stage 10j_plot_annotation_abundance_ec_after_soupx_harmony complete")
}

main()
