#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(scales)
  library(Cairo)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", paste0(..., collapse = ""), "\n")
}

main <- function() {
  stage_name <- "10ya_fix_overall_sample_contribution_stacked_colorlocked"
  input_rds <- "results/objects/08_unintegrated_clustered_umap.rds"
  sample_col <- "sample_id"
  output_dir <- "results/abundance"

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage ", stage_name)
  message_ts("Loaded object: ", input_rds)

  seu <- readRDS(input_rds)

  if (!sample_col %in% colnames(seu@meta.data)) {
    stop("Sample column not found: ", sample_col)
  }

  meta <- seu@meta.data

  # REQUIRED sample order everywhere in the presentation
  sample_order <- c("P7", "P12_CTRL", "P12_OIR", "P17_CTRL", "P17_OIR")

  # LOCKED sample-color mapping
  # If needed later, we only edit this vector and nothing else.
  sample_colors <- c(
    "P7"       = "#009E73",  # green
    "P12_CTRL" = "#F0E442",  # yellow
    "P12_OIR"  = "#0072B2",  # dark blue
    "P17_CTRL" = "#56B4E9",  # light blue
    "P17_OIR"  = "#E69F00"   # orange
  )

  observed_samples <- unique(as.character(meta[[sample_col]]))
  missing_samples <- setdiff(sample_order, observed_samples)
  extra_samples <- setdiff(observed_samples, sample_order)

  if (length(missing_samples) > 0) {
    stop("Missing expected sample IDs in object: ", paste(missing_samples, collapse = ", "))
  }
  if (length(extra_samples) > 0) {
    stop("Unexpected sample IDs in object: ", paste(extra_samples, collapse = ", "))
  }

  df <- meta %>%
    count(.data[[sample_col]], name = "n_cells") %>%
    mutate(
      sample_id = .data[[sample_col]],
      sample_id = factor(sample_id, levels = sample_order)
    ) %>%
    arrange(sample_id) %>%
    mutate(
      percent = 100 * n_cells / sum(n_cells),
      label = paste0(round(percent, 1), "%")
    )

  message_ts("Using sample order: ", paste(sample_order, collapse = ", "))
  message_ts("Using locked colors: ",
             paste(names(sample_colors), sample_colors, sep = "=", collapse = "; "))

  write_csv(
    df %>% select(sample_id, n_cells, percent),
    file.path(output_dir, "10ya_overall_sample_contribution_stacked_colorlocked_table.csv")
  )

  p <- ggplot(df, aes(x = "All samples", y = percent, fill = sample_id)) +
    geom_col(width = 0.6, color = NA) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5),
      size = 9,
      color = "black"
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      labels = function(x) paste0(x, "%"),
      expand = c(0, 0)
    ) +
    scale_fill_manual(
      values = sample_colors,
      breaks = sample_order,
      drop = FALSE,
      name = "Sample"
    ) +
    labs(
      title = "Overall cell contribution by sample",
      x = NULL,
      y = "Percent of total cells"
    ) +
    theme_classic(base_size = 20) +
    theme(
      plot.title = element_text(face = "bold", size = 28, hjust = 0.5),
      axis.title.y = element_text(face = "bold", size = 24),
      axis.text.x = element_text(size = 22),
      axis.text.y = element_text(size = 20, face = "bold"),
      axis.line = element_line(linewidth = 1.2, color = "black"),
      axis.ticks = element_line(linewidth = 1.2, color = "black"),
      legend.title = element_text(face = "bold", size = 18),
      legend.text = element_text(size = 16)
    )

  p_nolegend <- p + theme(legend.position = "none")

  CairoPDF(
    file.path(output_dir, "10ya_overall_sample_contribution_stacked_colorlocked.pdf"),
    width = 8,
    height = 10
  )
  print(p)
  dev.off()

  CairoPNG(
    file.path(output_dir, "10ya_overall_sample_contribution_stacked_colorlocked.png"),
    width = 8,
    height = 10,
    units = "in",
    dpi = 300
  )
  print(p)
  dev.off()

  CairoPDF(
    file.path(output_dir, "10ya_overall_sample_contribution_stacked_NOlegend_colorlocked.pdf"),
    width = 8,
    height = 10
  )
  print(p_nolegend)
  dev.off()

  CairoPNG(
    file.path(output_dir, "10ya_overall_sample_contribution_stacked_NOlegend_colorlocked.png"),
    width = 8,
    height = 10,
    units = "in",
    dpi = 300
  )
  print(p_nolegend)
  dev.off()

  message_ts("Stage ", stage_name, " complete")
}

main()
