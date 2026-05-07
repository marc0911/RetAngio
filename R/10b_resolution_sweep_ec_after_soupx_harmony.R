#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(clustree)
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

main <- function() {
  input_rds <- "results/objects/10a_ec_after_soupx_harmony_integrated.rds"
  output_rds <- "results/objects/10b_ec_after_soupx_harmony_resolution_sweep.rds"
  outdir <- "plots/10b_ec_after_soupx_harmony_resolution_sweep"

  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 10b_resolution_sweep_ec_after_soupx_harmony")

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  seu <- retangio_add_plotting_metadata(seu)

  if (!"harmony" %in% names(seu@reductions)) {
    stop("harmony reduction not found in object.")
  }

  resolutions <- seq(0.1, 1.0, by = 0.1)

  message_ts("Running FindNeighbors on harmony")
  seu <- FindNeighbors(
    object = seu,
    reduction = "harmony",
    dims = 1:30,
    verbose = TRUE
  )

  for (res in resolutions) {
    message_ts("Running FindClusters at resolution", res)
    seu <- FindClusters(
      object = seu,
      resolution = res,
      verbose = TRUE
    )
  }

  message_ts("Running UMAP on harmony")
  seu <- RunUMAP(
    object = seu,
    reduction = "harmony",
    dims = 1:30,
    reduction.name = "umap.harmony",
    reduction.key = "UMAPHARM_",
    verbose = TRUE
  )

  saveRDS(seu, output_rds)
  message_ts("Saved resolution sweep object:", output_rds)

  lims <- retangio_get_fixed_umap_limits(seu, "umap.harmony")

  if ("plot_condition" %in% colnames(seu@meta.data) && any(!is.na(seu$plot_condition))) {
    df_cond <- retangio_make_umap_df(
      obj = seu,
      reduction = "umap.harmony",
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
      title_text = "UMAP — Harmony — colored by condition",
      legend_title = "Condition",
      legend_ncol = 1
    )

    retangio_save_plot_pdf_png_safe(
      plot_obj = p_cond + theme(legend.position = "none"),
      filename_base = file.path(outdir, "10b_harmony_UMAP_byCondition_NOlegend"),
      width = retangio_width_main,
      height = retangio_height_main
    )

    p_cond_leg <- retangio_make_manual_legend_plot(
      labels = retangio_condition_labels[retangio_condition_levels],
      colors = setNames(retangio_condition_palette[retangio_condition_levels], retangio_condition_labels[retangio_condition_levels]),
      title_text = "Condition",
      ncol = 1
    )

    retangio_save_plot_pdf_png_safe(
      plot_obj = p_cond_leg,
      filename_base = file.path(outdir, "10b_harmony_UMAP_byCondition_legendOnly"),
      width = 4,
      height = 3
    )

    for (cond in retangio_condition_levels) {
      p_hi <- retangio_plot_umap_highlight(
        df = df_cond,
        group_col = "plot_condition",
        highlight_level = cond,
        palette_values = retangio_condition_palette,
        xlim = lims$xlim,
        ylim = lims$ylim,
        title_text = paste0("UMAP — Harmony — highlight ", retangio_condition_labels[[cond]]),
        legend_title = "Condition"
      ) + theme(legend.position = "none")

      retangio_save_plot_pdf_png_safe(
        plot_obj = p_hi,
        filename_base = file.path(outdir, paste0("10b_harmony_UMAP_highlight_", cond)),
        width = retangio_width_main,
        height = retangio_height_main
      )
    }
  }

  for (res in resolutions) {
    cluster_col <- paste0("SCT_snn_res.", res)
    cluster_levels_num <- sort(unique(as.numeric(as.character(seu@meta.data[[cluster_col]]))))
    cluster_levels <- as.character(cluster_levels_num)
    cluster_pal <- retangio_make_cluster_palette(cluster_levels)

    df_cluster <- retangio_make_umap_df(
      obj = seu,
      reduction = "umap.harmony",
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
      title_text = paste0("UMAP — Harmony — colored by cluster — res ", sprintf("%.1f", res)),
      legend_title = "Cluster",
      legend_ncol = 2
    )

    stem <- file.path(outdir, paste0("10b_harmony_res_", sprintf("%.1f", res), "_UMAP_byCluster"))

    retangio_save_plot_pdf_png_safe(
      plot_obj = p_cluster + theme(legend.position = "none"),
      filename_base = paste0(stem, "_NOlegend"),
      width = retangio_width_main,
      height = retangio_height_main
    )

    p_cluster_leg <- retangio_make_manual_legend_plot(
      labels = cluster_levels,
      colors = cluster_pal,
      title_text = "Cluster",
      ncol = 2
    )

    legend_height <- max(4.5, ceiling(ceiling(length(cluster_levels) / 2)) * 0.45 + 1.6)

    retangio_save_plot_pdf_png_safe(
      plot_obj = p_cluster_leg,
      filename_base = paste0(stem, "_legendOnly"),
      width = 5.8,
      height = legend_height
    )
  }

  message_ts("Generating clustree")
  p_tree <- clustree::clustree(
    seu@meta.data,
    prefix = "SCT_snn_res."
  ) +
    theme_classic(base_size = 16, base_family = "Arial") +
    labs(title = "Clustree — Harmony EC after SoupX")

  retangio_save_plot_pdf_png_safe(
    plot_obj = p_tree,
    filename_base = file.path(outdir, "10b_harmony_clustree"),
    width = 12,
    height = 10
  )

  message_ts("Stage 10b_resolution_sweep_ec_after_soupx_harmony complete")
}

main()
