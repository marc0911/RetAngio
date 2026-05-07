# ------------------------------------------------------------------
# RetAngio presentation-only plotting helper
# ------------------------------------------------------------------
# IMPORTANT:
# This helper was created for late presentation/export scripts only.
# It is currently used by:
#   - scripts/10x_presentation_plots_tipcells_ctrl_only.R
#   - scripts/10x_presentation_plots_tipcells_overall.R
#
# Canonical helper for the main analysis pipeline:
#   - R/plotting_helpers.R
#
# Do not merge or delete this file until the presentation scripts are
# either retired or refactored in a controlled cleanup step.
# ------------------------------------------------------------------

#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(cowplot)
  library(Cairo)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

retangio_pt_size <- 0.10
retangio_pt_alpha <- 0.45
retangio_base_font <- 20

retangio_condition_levels <- c("P7", "P12_CTRL", "P12_OIR", "P17_CTRL", "P17_OIR")
retangio_condition_labels <- c(
  "P7" = "P7",
  "P12_CTRL" = "P12 Ctrl",
  "P12_OIR" = "P12 OIR",
  "P17_CTRL" = "P17 Ctrl",
  "P17_OIR" = "P17 OIR"
)

retangio_condition_palette <- c(
  "P7" = "#E69F00",
  "P12_CTRL" = "#56B4E9",
  "P12_OIR" = "#009E73",
  "P17_CTRL" = "#0072B2",
  "P17_OIR" = "#D55E00"
)

retangio_theme_proj <- function() {
  theme_classic(base_size = retangio_base_font, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 16),
      legend.key.height = unit(0.9, "cm"),
      legend.key.width = unit(0.9, "cm")
    )
}

retangio_save_plot <- function(plot_obj, filename_base, width, height, dpi = 300) {
  Cairo::CairoPDF(
    file = paste0(filename_base, ".pdf"),
    width = width,
    height = height,
    family = "Arial"
  )
  print(plot_obj)
  dev.off()

  Cairo::CairoPNG(
    filename = paste0(filename_base, ".png"),
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  print(plot_obj)
  dev.off()
}

retangio_extract_legend_plot <- function(plot_obj) {
  leg <- cowplot::get_legend(plot_obj + theme(legend.position = "right"))
  cowplot::ggdraw(leg)
}

retangio_guess_plot_condition <- function(meta_df) {
  if ("plot_condition" %in% colnames(meta_df)) {
    vals <- as.character(meta_df$plot_condition)
    if (sum(!is.na(vals)) > 0) {
      return(factor(vals, levels = retangio_condition_levels))
    }
  }

  if ("sample_id" %in% colnames(meta_df)) {
    vals <- as.character(meta_df$sample_id)
    return(factor(vals, levels = retangio_condition_levels))
  }

  if ("orig.ident" %in% colnames(meta_df)) {
    vals <- as.character(meta_df$orig.ident)
    return(factor(vals, levels = retangio_condition_levels))
  }

  stop("Could not determine plot_condition from metadata.")
}

retangio_add_plotting_metadata <- function(obj) {
  obj$plot_condition <- retangio_guess_plot_condition(obj@meta.data)
  obj
}

retangio_detect_reduction <- function(obj, preferred = character(0)) {
  reds <- Reductions(obj)

  for (r in preferred) {
    if (r %in% reds) {
      return(r)
    }
  }

  umaps <- grep("^umap", reds, value = TRUE)
  if (length(umaps) > 0) {
    return(umaps[1])
  }

  stop(
    "No suitable UMAP reduction found. Available reductions: ",
    paste(reds, collapse = ", ")
  )
}

retangio_get_fixed_umap_limits <- function(obj, reduction) {
  emb <- Embeddings(obj, reduction = reduction)
  list(
    xlim = range(emb[, 1]),
    ylim = range(emb[, 2])
  )
}

retangio_make_umap_df <- function(obj, reduction, color_col, levels_use = NULL, shuffle_seed = 1) {
  emb <- Embeddings(obj, reduction = reduction)
  df <- as.data.frame(emb)
  colnames(df)[1:2] <- c("UMAP_1", "UMAP_2")
  df[[color_col]] <- obj@meta.data[rownames(df), color_col, drop = TRUE]

  if (!is.null(levels_use)) {
    df[[color_col]] <- factor(as.character(df[[color_col]]), levels = levels_use)
  }

  set.seed(shuffle_seed)
  df <- df[sample(nrow(df)), , drop = FALSE]
  df
}

retangio_plot_umap_discrete <- function(
  df,
  color_col,
  palette_values,
  xlim,
  ylim,
  title_text,
  legend_title
) {
  ggplot(df, aes(UMAP_1, UMAP_2, color = .data[[color_col]])) +
    geom_point(size = retangio_pt_size, alpha = retangio_pt_alpha) +
    scale_color_manual(
      values = palette_values,
      breaks = levels(df[[color_col]]),
      limits = levels(df[[color_col]]),
      drop = FALSE,
      name = legend_title
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = title_text, x = "UMAP 1", y = "UMAP 2") +
    retangio_theme_proj() +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))
}

retangio_sort_cluster_levels <- function(x) {
  x <- unique(as.character(x))
  x_num <- suppressWarnings(as.numeric(x))

  if (all(!is.na(x_num))) {
    x[order(x_num)]
  } else {
    sort(x)
  }
}

retangio_make_cluster_palette <- function(cluster_levels) {
  pal <- scales::hue_pal(h = c(0, 360) + 15, c = 100, l = 65)(length(cluster_levels))
  names(pal) <- cluster_levels
  pal
}

retangio_make_abundance_df <- function(obj, cluster_col) {
  if (!cluster_col %in% colnames(obj@meta.data)) {
    stop("Cluster column not found: ", cluster_col)
  }

  meta <- obj@meta.data
  meta$plot_condition <- obj$plot_condition

  df <- meta %>%
    dplyr::transmute(
      plot_condition = as.character(plot_condition),
      cluster = as.character(.data[[cluster_col]])
    ) %>%
    dplyr::filter(!is.na(plot_condition), !is.na(cluster)) %>%
    dplyr::group_by(plot_condition, cluster) %>%
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop")

  sample_levels_present <- retangio_condition_levels[retangio_condition_levels %in% unique(df$plot_condition)]
  cluster_levels <- retangio_sort_cluster_levels(df$cluster)

  df <- df %>%
    dplyr::mutate(
      plot_condition = factor(plot_condition, levels = sample_levels_present),
      cluster = factor(cluster, levels = cluster_levels)
    ) %>%
    dplyr::group_by(plot_condition) %>%
    dplyr::mutate(
      percent = n_cells / sum(n_cells),
      pct_label = ifelse(percent >= 0.05, scales::percent(percent, accuracy = 0.1), "")
    ) %>%
    dplyr::ungroup()

  df
}

retangio_plot_abundance <- function(df, title_text, legend_title) {
  cluster_levels <- levels(df$cluster)
  cluster_pal <- retangio_make_cluster_palette(cluster_levels)

  x_labels <- retangio_condition_labels[levels(df$plot_condition)]

  ggplot(df, aes(x = plot_condition, y = percent, fill = cluster)) +
    geom_col(width = 0.82, color = NA) +
    geom_text(
      aes(label = pct_label),
      position = position_stack(vjust = 0.5),
      size = 5,
      color = "black"
    ) +
    scale_fill_manual(
      values = cluster_pal,
      breaks = cluster_levels,
      limits = cluster_levels,
      drop = FALSE,
      name = legend_title
    ) +
    scale_x_discrete(labels = x_labels) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    labs(
      title = title_text,
      x = NULL,
      y = "Percent of sample cells"
    ) +
    retangio_theme_proj() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    guides(fill = guide_legend(override.aes = list(size = 6, alpha = 1)))
}
