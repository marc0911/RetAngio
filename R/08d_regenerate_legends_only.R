#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

condition_levels <- c("P7", "P12_CTRL", "P12_OIR", "P17_CTRL", "P17_OIR")
condition_labels <- c("P7", "P12 Ctrl", "P12 OIR", "P17 Ctrl", "P17 OIR")
names(condition_labels) <- condition_levels

condition_palette <- c(
  "P7" = "#E69F00",
  "P12_CTRL" = "#56B4E9",
  "P12_OIR" = "#009E73",
  "P17_CTRL" = "#0072B2",
  "P17_OIR" = "#D55E00"
)

make_cluster_palette <- function(cluster_levels) {
  cluster_levels <- as.character(cluster_levels)
  pal <- scales::hue_pal(h = c(0, 360) + 15, c = 100, l = 65)(length(cluster_levels))
  names(pal) <- cluster_levels
  pal
}

save_plot_pdf_png <- function(plot_obj, stem, width, height, dpi = 300) {
  grDevices::cairo_pdf(
    filename = paste0(stem, ".pdf"),
    width = width,
    height = height,
    family = "Arial"
  )
  print(plot_obj)
  dev.off()

  png(
    filename = paste0(stem, ".png"),
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

make_manual_legend_plot <- function(labels, colors, title_text, ncol = 1) {
  labels <- as.character(labels)

  if (ncol == 1) {
    df <- data.frame(
      label = labels,
      color = unname(colors[labels]),
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

  # ncol = 2 for cluster legends
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
  input_rds <- "results/objects/08d_ec_after_soupx_unintegrated_resolution_sweep.rds"
  outdir <- "plots/08d_ec_after_soupx_unintegrated_resolution_sweep"

  message_ts("Starting legends-only regeneration")
  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  # Condition legend
  if ("plot_condition" %in% colnames(seu@meta.data)) {
    present_conditions <- intersect(condition_levels, unique(as.character(seu@meta.data$plot_condition)))
    if (length(present_conditions) > 0) {
      p_cond <- make_manual_legend_plot(
        labels = condition_labels[present_conditions],
        colors = setNames(condition_palette[present_conditions], condition_labels[present_conditions]),
        title_text = "Condition",
        ncol = 1
      )

      save_plot_pdf_png(
        plot_obj = p_cond,
        stem = file.path(outdir, "08d_unintegrated_UMAP_byCondition_legendOnly"),
        width = 4,
        height = 3
      )
      message_ts("Condition legend regenerated")
    }
  }

  # Cluster legends for all resolutions
  res_cols <- grep("^SCT_snn_res\\.", colnames(seu@meta.data), value = TRUE)
  res_cols <- sort(res_cols)

  for (cluster_col in res_cols) {
    res_value <- sub("^SCT_snn_res\\.", "", cluster_col)
    cluster_vals <- as.character(seu@meta.data[[cluster_col]])
    cluster_levels <- as.character(sort(unique(as.numeric(cluster_vals))))
    cluster_palette <- make_cluster_palette(cluster_levels)

    p_cluster <- make_manual_legend_plot(
      labels = cluster_levels,
      colors = cluster_palette,
      title_text = "Cluster",
      ncol = 2
    )

    stem <- file.path(
      outdir,
      paste0("08d_unintegrated_res_", sprintf("%.1f", as.numeric(res_value)), "_UMAP_byCluster_legendOnly")
    )

    legend_height <- max(4.5, ceiling(ceiling(length(cluster_levels) / 2)) * 0.45 + 1.6)

    save_plot_pdf_png(
      plot_obj = p_cluster,
      stem = stem,
      width = 5.8,
      height = legend_height
    )

    message_ts("Cluster legend regenerated for resolution", res_value)
  }

  message_ts("Legends-only regeneration complete")
}

main()
