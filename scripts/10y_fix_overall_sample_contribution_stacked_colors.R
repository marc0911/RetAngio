suppressPackageStartupMessages({
  .libPaths(c("/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs", .libPaths()))
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(readr)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

choose_sample_col <- function(meta) {
  candidates <- c("plot_condition", "sample_id", "orig.ident")
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0) {
    stop(
      "Could not find a sample column. Available columns: ",
      paste(colnames(meta), collapse = ", ")
    )
  }
  hit[1]
}

main <- function() {
  message_ts("Starting stage 10y_fix_overall_sample_contribution_stacked_colors")

  input_rds <- "results/objects/08_unintegrated_clustered_umap.rds"
  output_dir <- "results/abundance"

  if (!file.exists(input_rds)) {
    stop("Input object not found: ", input_rds)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  seu <- readRDS(input_rds)
  message_ts("Loaded object: ", input_rds)

  meta <- seu@meta.data
  sample_col <- choose_sample_col(meta)
  message_ts("Using sample column: ", sample_col)

  sample_levels <- c("P7", "P12_CTRL", "P12_OIR", "P17_CTRL", "P17_OIR")
  sample_colors <- c(
    "P7" = "#E69F00",
    "P12_CTRL" = "#56B4E9",
    "P12_OIR" = "#0072B2",
    "P17_CTRL" = "#F0E442",
    "P17_OIR" = "#009E73"
  )

  counts_df <- meta %>%
    mutate(sample_plot = .data[[sample_col]]) %>%
    count(sample_plot, name = "n_cells") %>%
    mutate(sample_plot = as.character(sample_plot)) %>%
    filter(sample_plot %in% sample_levels) %>%
    mutate(
      sample_plot = factor(sample_plot, levels = sample_levels),
      percent = n_cells / sum(n_cells),
      label = sprintf("%.1f%%", percent * 100)
    ) %>%
    arrange(sample_plot)

  write_csv(
    counts_df %>% mutate(sample_plot = as.character(sample_plot)),
    file.path(output_dir, "10y_overall_sample_contribution_stacked_corrected_table.csv")
  )

  p <- ggplot(
    counts_df,
    aes(x = "All samples", y = percent, fill = sample_plot)
  ) +
    geom_col(
      width = 0.55,
      position = position_stack(reverse = TRUE),
      color = NA
    ) +
    geom_text(
      aes(label = label),
      position = position_stack(vjust = 0.5, reverse = TRUE),
      size = 10
    ) +
    scale_fill_manual(values = sample_colors, breaks = sample_levels, drop = FALSE) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    labs(
      title = "Overall cell contribution by sample",
      x = NULL,
      y = "Percent of total cells"
    ) +
    theme_classic(base_size = 28) +
    theme(
      legend.position = "none",
      axis.title.y = element_text(face = "bold", size = 28),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 24),
      axis.text.y = element_text(size = 24),
      axis.line = element_line(linewidth = 1.5),
      axis.ticks = element_line(linewidth = 1.5),
      axis.ticks.length = unit(0.25, "cm"),
      plot.title = element_text(face = "bold", size = 34, hjust = 0.5),
      plot.margin = margin(15, 15, 15, 15)
    )

  pdf(
    file.path(output_dir, "10y_overall_sample_contribution_stacked_NOlegend_corrected.pdf"),
    width = 10,
    height = 10
  )
  print(p)
  dev.off()

  png(
    filename = file.path(output_dir, "10y_overall_sample_contribution_stacked_NOlegend_corrected.png"),
    width = 3000,
    height = 3000,
    res = 300,
    type = "cairo"
  )
  print(p)
  dev.off()

  message_ts("Stage 10y_fix_overall_sample_contribution_stacked_colors complete")
}

main()
