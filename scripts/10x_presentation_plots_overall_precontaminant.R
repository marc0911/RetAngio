#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(cowplot)
  library(Cairo)
})

source("R/plotting_helpers.R")

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

retangio_save_plot_png_base <- function(plot_obj, filename_base, width, height, dpi = 300) {
  png(
    filename = paste0(filename_base, ".png"),
    width = width,
    height = height,
    units = "in",
    res = dpi,
    bg = "white",
    type = "cairo"
  )
  print(plot_obj)
  dev.off()
}

retangio_save_plot_pdf_png_safe <- function(plot_obj, filename_base, width, height, dpi = 300) {
  grDevices::cairo_pdf(
    filename = paste0(filename_base, ".pdf"),
    width = width,
    height = height,
    family = "Arial"
  )
  print(plot_obj)
  dev.off()

  retangio_save_plot_png_base(
    plot_obj = plot_obj,
    filename_base = filename_base,
    width = width,
    height = height,
    dpi = dpi
  )
}

retangio_make_manual_legend_plot <- function(labels, colors, title_text, ncol = 1) {
  labels <- as.character(labels)

  if (ncol == 1) {
    df <- data.frame(
      label = labels,
      x = 1,
      y = rev(seq_along(labels)),
      stringsAsFactors = FALSE
    )

    p <- ggplot(df, aes(x = x, y = y)) +
      geom_point(aes(color = label), size = 5) +
      geom_text(aes(label = label), hjust = 0, nudge_x = 0.25, family = "Arial", size = 5) +
      scale_color_manual(values = colors[labels], breaks = labels, limits = labels, drop = FALSE) +
      labs(title = title_text) +
      xlim(0.7, 3.5) +
      ylim(0.5, length(labels) + 0.5) +
      theme_void(base_size = 18) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0, family = "Arial", size = 18),
        legend.position = "none",
        plot.margin = margin(15, 15, 15, 15)
      )

    return(p)
  }

  n <- length(labels)
  n_left <- ceiling(n / 2)
  left_labels <- labels[1:n_left]
  right_labels <- labels[(n_left + 1):n]

  df_left <- data.frame(
    label = left_labels,
    x = 1,
    y = rev(seq_along(left_labels)),
    stringsAsFactors = FALSE
  )

  df_right <- NULL
  if (length(right_labels) > 0) {
    df_right <- data.frame(
      label = right_labels,
      x = 4,
      y = rev(seq_along(right_labels)),
      stringsAsFactors = FALSE
    )
  }

  df <- rbind(df_left, df_right)
  max_y <- max(df$y)

  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(aes(color = label), size = 5) +
    geom_text(aes(label = label), hjust = 0, nudge_x = 0.25, family = "Arial", size = 5) +
    scale_color_manual(values = colors[labels], breaks = labels, limits = labels, drop = FALSE) +
    labs(title = title_text) +
    xlim(0.7, 6.5) +
    ylim(0.5, max_y + 0.5) +
    theme_void(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, family = "Arial", size = 18),
      legend.position = "none",
      plot.margin = margin(15, 15, 15, 15)
    )

  p
}

detect_reduction <- function(seu) {
  candidates <- c("umap.unintegrated", "umap.integrated", "umap.rpca", "umap.harmony", "umap")
  hit <- candidates[candidates %in% names(seu@reductions)]
  if (length(hit) == 0) {
    stop("No suitable UMAP reduction found. Available reductions: ",
         paste(names(seu@reductions), collapse = ", "))
  }
  hit[1]
}

detect_cluster_col <- function(seu) {
  if ("seurat_clusters" %in% colnames(seu@meta.data)) return("seurat_clusters")

  cluster_cols <- grep("^SCT_snn_res\\.", colnames(seu@meta.data), value = TRUE)
  if (length(cluster_cols) > 0) return(tail(sort(cluster_cols), 1))

  ids <- as.character(Idents(seu))
  if (length(ids) == ncol(seu)) {
    seu$seurat_clusters <- ids
    return("seurat_clusters")
  }

  stop("Could not detect a cluster column.")
}

sort_cluster_levels <- function(x) {
  x_chr <- as.character(x)
  x_num <- suppressWarnings(as.numeric(x_chr))
  if (all(!is.na(x_num))) return(as.character(sort(unique(x_num))))
  sort(unique(x_chr))
}

main <- function() {
  input_rds <- "results/objects/08_unintegrated_clustered_umap.rds"
  outdir <- "plots/10x_presentation_overall_precontaminant"
  abundance_csv <- "results/abundance/10x_presentation_overall_precontaminant_sample_contribution.csv"

  message_ts("Starting stage 10x_presentation_plots_overall_precontaminant")

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  seu <- retangio_add_plotting_metadata(seu)

  reduction_use <- detect_reduction(seu)
  cluster_col <- detect_cluster_col(seu)

  message_ts("Using reduction:", reduction_use)
  message_ts("Using cluster column:", cluster_col)

  if (!"plot_condition" %in% colnames(seu@meta.data) || all(is.na(seu$plot_condition))) {
    stop("plot_condition could not be inferred from sample_id/orig.ident.")
  }

  lims <- retangio_get_fixed_umap_limits(seu, reduction_use)

  cluster_levels <- sort_cluster_levels(seu@meta.data[[cluster_col]])
  cluster_pal <- retangio_make_cluster_palette(cluster_levels)

  df_cluster <- retangio_make_umap_df(
    obj = seu,
    reduction = reduction_use,
    color_col = cluster_col,
    levels_use = cluster_levels,
    shuffle_seed = 1
  )

  p_cluster <- retangio_plot_umap_discrete(
    df = df_cluster,
    color_col = cluster_col,
    palette_values = cluster_pal,
    xlim = lims$xlim,
    ylim = lims$ylim,
    title_text = "UMAP — overall object — colored by cluster",
    legend_title = "Cluster",
    legend_ncol = 2
  )

  retangio_save_plot_pdf_png_safe(
    plot_obj = p_cluster + theme(legend.position = "none"),
    filename_base = file.path(outdir, "10x_overall_UMAP_byCluster_NOlegend"),
    width = retangio_width_main,
    height = retangio_height_main
  )

  p_cluster_leg <- retangio_make_manual_legend_plot(
    labels = cluster_levels,
    colors = cluster_pal,
    title_text = "Cluster",
    ncol = 2
  )

  legend_height_cluster <- max(4.5, ceiling(ceiling(length(cluster_levels) / 2)) * 0.45 + 1.6)

  retangio_save_plot_pdf_png_safe(
    plot_obj = p_cluster_leg,
    filename_base = file.path(outdir, "10x_overall_UMAP_byCluster_legendOnly"),
    width = 5.8,
    height = legend_height_cluster
  )

  df_cond <- retangio_make_umap_df(
    obj = seu,
    reduction = reduction_use,
    color_col = "plot_condition",
    levels_use = retangio_condition_levels,
    shuffle_seed = 1
  )

  p_cond <- retangio_plot_umap_discrete(
    df = df_cond,
    color_col = "plot_condition",
    palette_values = retangio_condition_palette,
    xlim = lims$xlim,
    ylim = lims$ylim,
    title_text = "UMAP — overall object — colored by sample",
    legend_title = "Sample",
    legend_ncol = 1
  )

  retangio_save_plot_pdf_png_safe(
    plot_obj = p_cond + theme(legend.position = "none"),
    filename_base = file.path(outdir, "10x_overall_UMAP_bySample_NOlegend"),
    width = retangio_width_main,
    height = retangio_height_main
  )

  p_cond_leg <- retangio_make_manual_legend_plot(
    labels = retangio_condition_labels[retangio_condition_levels],
    colors = setNames(
      retangio_condition_palette[retangio_condition_levels],
      retangio_condition_labels[retangio_condition_levels]
    ),
    title_text = "Sample",
    ncol = 1
  )

  retangio_save_plot_pdf_png_safe(
    plot_obj = p_cond_leg,
    filename_base = file.path(outdir, "10x_overall_UMAP_bySample_legendOnly"),
    width = 4,
    height = 3
  )

  abundance_df <- seu@meta.data %>%
    mutate(plot_condition = factor(as.character(plot_condition), levels = retangio_condition_levels)) %>%
    count(plot_condition, name = "n_cells", .drop = FALSE) %>%
    mutate(
      fraction = n_cells / sum(n_cells),
      percent = 100 * fraction,
      sample_label = retangio_condition_labels[as.character(plot_condition)],
      percent_label = scales::percent(fraction, accuracy = 0.1)
    )

  write.csv(abundance_df, abundance_csv, row.names = FALSE)

  ymax <- max(abundance_df$fraction)
  ymax <- max(0.05, ymax * 1.18)

  p_abund <- ggplot(abundance_df, aes(x = plot_condition, y = fraction, fill = plot_condition)) +
    geom_col(width = 0.8) +
    geom_text(
      aes(label = percent_label),
      vjust = -0.35,
      family = "Arial",
      size = 5
    ) +
    scale_fill_manual(
      values = retangio_condition_palette,
      breaks = retangio_condition_levels,
      limits = retangio_condition_levels,
      drop = FALSE
    ) +
    scale_x_discrete(
      limits = retangio_condition_levels,
      labels = retangio_condition_labels[retangio_condition_levels],
      drop = FALSE
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, ymax),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = "Overall cells — contribution of each sample",
      x = NULL,
      y = "Fraction of all cells"
    ) +
    retangio_theme_proj() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  retangio_save_plot_pdf_png_safe(
    plot_obj = p_abund,
    filename_base = file.path(outdir, "10x_overall_sample_contribution_NOlegend"),
    width = retangio_width_main,
    height = retangio_height_main
  )

  message_ts("Saved abundance table:", abundance_csv)
  message_ts("Stage 10x_presentation_plots_overall_precontaminant complete")
}

main()
