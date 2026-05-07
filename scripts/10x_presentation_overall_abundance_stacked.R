#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

log_message <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_rds <- "results/objects/08_unintegrated_clustered_umap.rds"
  out_plot_dir <- "plots/10x_presentation_overall_precontaminant"
  out_table <- "results/abundance/10x_presentation_overall_precontaminant_sample_contribution.csv"

  dir.create(out_plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(out_table), recursive = TRUE, showWarnings = FALSE)

  log_message("Starting stage 10x_presentation_overall_abundance_stacked")

  seu <- readRDS(input_rds)
  log_message("Loaded object:", input_rds)

  meta <- seu@meta.data

  sample_col <- NULL
  if ("sample_id" %in% colnames(meta)) {
    sample_col <- "sample_id"
  } else if ("orig.ident" %in% colnames(meta)) {
    sample_col <- "orig.ident"
  } else {
    stop("No sample column found (expected sample_id or orig.ident).")
  }

  sample_order <- c("P7", "P12_CTRL", "P12_OIR", "P17_CTRL", "P17_OIR")

  ## Set colors here if you want to tweak them later
  sample_colors <- c(
    "P7"       = "#E69F00",
    "P12_CTRL" = "#56B4E9",
    "P12_OIR"  = "#0072B2",
    "P17_CTRL" = "#F0E442",
    "P17_OIR"  = "#009E73"
  )

  abundance_df <- meta |>
    dplyr::count(.data[[sample_col]], name = "n_cells") |>
    dplyr::rename(sample = 1)

  abundance_df$sample <- factor(abundance_df$sample, levels = sample_order)

  abundance_df <- abundance_df |>
    dplyr::arrange(sample) |>
    dplyr::mutate(
      percent = 100 * n_cells / sum(n_cells),
      label = sprintf("%.1f%%", percent),
      bar = "All samples"
    )

  write.csv(abundance_df, out_table, row.names = FALSE)

  p <- ggplot(abundance_df, aes(x = bar, y = percent, fill = sample)) +
    geom_col(width = 0.62, color = NA) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 8
    ) +
    scale_fill_manual(
      values = sample_colors,
      breaks = sample_order,
      drop = FALSE
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      labels = function(x) paste0(x, "%"),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = "Overall cell contribution by sample",
      x = NULL,
      y = "Percent of total cells"
    ) +
    theme_classic(base_size = 22) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 30, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 24),
      axis.text = element_text(size = 20, colour = "black"),
      axis.line = element_line(linewidth = 1.4, colour = "black"),
      axis.ticks = element_line(linewidth = 1.4, colour = "black")
    )

  ggsave(
    filename = file.path(out_plot_dir, "10x_presentation_overall_sample_contribution_stacked_NOlegend.png"),
    plot = p,
    width = 10,
    height = 8,
    dpi = 300
  )

  ggsave(
    filename = file.path(out_plot_dir, "10x_presentation_overall_sample_contribution_stacked_NOlegend.pdf"),
    plot = p,
    width = 10,
    height = 8
  )

  log_message("Stage 10x_presentation_overall_abundance_stacked complete")
}

main()
