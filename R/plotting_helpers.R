#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(cowplot)
  library(Cairo)
})

# ------------------------------------------------------------------
# Global plotting conventions
# ------------------------------------------------------------------

retangio_pt_size <- 0.10
retangio_pt_alpha <- 0.45
retangio_bg_col <- "grey80"
retangio_bg_alpha <- 0.12
retangio_bg_size <- 0.10
retangio_base_font <- 20
retangio_width_main <- 11
retangio_height_main <- 8.5

retangio_condition_levels <- c("P7", "P12_CTRL", "P12_OIR", "P17_CTRL", "P17_OIR")
retangio_condition_labels <- c("P7", "P12 Ctrl", "P12 OIR", "P17 Ctrl", "P17 OIR")
names(retangio_condition_labels) <- retangio_condition_levels

retangio_condition_palette <- c(
  "P7" = "#E69F00",
  "P12_CTRL" = "#56B4E9",
  "P12_OIR" = "#009E73",
  "P17_CTRL" = "#0072B2",
  "P17_OIR" = "#D55E00"
)

retangio_okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000"
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


retangio_save_plot_pdf_only <- function(plot_obj, filename_base, width, height) {
  Cairo::CairoPDF(
    file = paste0(filename_base, ".pdf"),
    width = width,
    height = height,
    family = "Arial"
  )
  print(plot_obj)
  dev.off()
}

retangio_make_named_palette <- function(levels_vec, preferred = NULL) {
  levels_vec <- as.character(levels_vec)

  if (is.null(preferred)) {
    preferred <- character(0)
  } else {
    preferred <- preferred[names(preferred) %in% levels_vec]
  }

  missing_levels <- setdiff(levels_vec, names(preferred))
  if (length(missing_levels) > 0) {
    extra_cols <- scales::hue_pal(h = c(0, 360) + 15, c = 100, l = 65)(length(missing_levels))
    names(extra_cols) <- missing_levels
    preferred <- c(preferred, extra_cols)
  }

  preferred[levels_vec]
}

retangio_make_cluster_palette <- function(cluster_levels) {
  cluster_levels <- as.character(cluster_levels)
  pal <- scales::hue_pal(h = c(0, 360) + 15, c = 100, l = 65)(length(cluster_levels))
  names(pal) <- cluster_levels
  pal
}

retangio_guess_plot_condition <- function(meta_df) {
  n <- nrow(meta_df)

  if ("sample_id" %in% colnames(meta_df)) {
    vals <- as.character(meta_df$sample_id)
    if (length(vals) == n && all(unique(vals) %in% retangio_condition_levels)) {
      return(factor(vals, levels = retangio_condition_levels))
    }
  }

  if ("orig.ident" %in% colnames(meta_df)) {
    vals <- as.character(meta_df$orig.ident)
    if (length(vals) == n && all(unique(vals) %in% retangio_condition_levels)) {
      return(factor(vals, levels = retangio_condition_levels))
    }
  }

  rep(NA, n)
}

retangio_add_plotting_metadata <- function(obj) {
  obj$plot_condition <- retangio_guess_plot_condition(obj@meta.data)
  obj
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
  legend_title,
  legend_ncol = 1
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
    guides(color = guide_legend(
      override.aes = list(size = 5, alpha = 1),
      ncol = legend_ncol
    ))
}

retangio_plot_umap_highlight <- function(
  df,
  group_col,
  highlight_level,
  palette_values,
  xlim,
  ylim,
  title_text,
  legend_title
) {
  df_hi <- df %>%
    mutate(.is_hi = (.data[[group_col]] == highlight_level))

  ggplot(df_hi, aes(UMAP_1, UMAP_2)) +
    geom_point(color = retangio_bg_col, size = retangio_bg_size, alpha = retangio_bg_alpha) +
    geom_point(
      data = df_hi %>% filter(.is_hi),
      aes(color = .data[[group_col]]),
      size = retangio_pt_size,
      alpha = 0.95
    ) +
    scale_color_manual(
      values = palette_values,
      breaks = levels(df[[group_col]]),
      limits = levels(df[[group_col]]),
      drop = FALSE,
      name = legend_title
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = title_text, x = "UMAP 1", y = "UMAP 2") +
    retangio_theme_proj()
}

retangio_extract_legend <- function(p) {
  cowplot::get_legend(p + theme(legend.position = "right"))
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

retangio_make_manual_legend_plot <- function(labels, colors, title_text, ncol = 1) {
  labels <- as.character(labels)

  if (ncol == 1) {
    df <- data.frame(
      label = labels,
      x = 1,
      y = rev(seq_along(labels)),
      stringsAsFactors = FALSE
    )

    return(
      ggplot(df, aes(x = x, y = y)) +
        geom_point(aes(color = label), size = 5) +
        geom_text(aes(label = label), hjust = 0, nudge_x = 0.25, family = "Arial", size = 5) +
        scale_color_manual(values = colors[labels], breaks = labels, limits = labels, drop = FALSE) +
        labs(title = title_text) +
        xlim(0.7, 3.8) +
        ylim(0.5, length(labels) + 0.5) +
        theme_void(base_size = 18) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0, family = "Arial", size = 18),
          legend.position = "none",
          plot.margin = margin(15, 15, 15, 15)
        )
    )
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

  ggplot(df, aes(x = x, y = y)) +
    geom_point(aes(color = label), size = 5) +
    geom_text(aes(label = label), hjust = 0, nudge_x = 0.25, family = "Arial", size = 5) +
    scale_color_manual(values = colors[labels], breaks = labels, limits = labels, drop = FALSE) +
    labs(title = title_text) +
    xlim(0.7, 6.6) +
    ylim(0.5, max_y + 0.5) +
    theme_void(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, family = "Arial", size = 18),
      legend.position = "none",
      plot.margin = margin(15, 15, 15, 15)
    )
}

retangio_plot_umap_highlight_subset_by_annotation <- function(
  obj,
  reduction,
  subset_col,
  subset_level,
  color_col,
  color_levels,
  palette_values,
  xlim,
  ylim,
  title_text
) {
  emb <- Embeddings(obj, reduction = reduction)
  df <- as.data.frame(emb)
  colnames(df)[1:2] <- c("UMAP_1", "UMAP_2")

  df[[subset_col]] <- obj@meta.data[rownames(df), subset_col, drop = TRUE]
  df[[color_col]] <- obj@meta.data[rownames(df), color_col, drop = TRUE]
  df[[color_col]] <- factor(as.character(df[[color_col]]), levels = color_levels)
  df$.is_hi <- as.character(df[[subset_col]]) == subset_level

  set.seed(1)
  df <- df[sample(nrow(df)), , drop = FALSE]

  ggplot(df, aes(UMAP_1, UMAP_2)) +
    geom_point(color = retangio_bg_col, size = retangio_bg_size, alpha = retangio_bg_alpha) +
    geom_point(
      data = df[df$.is_hi, , drop = FALSE],
      aes(color = .data[[color_col]]),
      size = retangio_pt_size,
      alpha = 0.95
    ) +
    scale_color_manual(
      values = palette_values,
      breaks = color_levels,
      limits = color_levels,
      drop = FALSE,
      name = "Annotation"
    ) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(title = title_text, x = "UMAP 1", y = "UMAP 2") +
    retangio_theme_proj()
}

retangio_plot_annotation_bundle <- function(
  obj,
  reduction,
  color_col,
  color_levels,
  palette_values,
  outdir,
  stem,
  title_text,
  legend_title = "Annotation",
  split_col = "plot_condition",
  split_levels = retangio_condition_levels,
  make_split = TRUE,
  make_highlight = TRUE
) {
  obj@meta.data[[color_col]] <- factor(as.character(obj@meta.data[[color_col]]), levels = color_levels)
  if (split_col %in% colnames(obj@meta.data)) {
    obj@meta.data[[split_col]] <- factor(as.character(obj@meta.data[[split_col]]), levels = split_levels)
  }

  lims <- retangio_get_fixed_umap_limits(obj, reduction)

  df <- retangio_make_umap_df(
    obj = obj,
    reduction = reduction,
    color_col = color_col,
    levels_use = color_levels,
    shuffle_seed = 1
  )

  p <- retangio_plot_umap_discrete(
    df = df,
    color_col = color_col,
    palette_values = palette_values,
    xlim = lims$xlim,
    ylim = lims$ylim,
    title_text = title_text,
    legend_title = legend_title,
    legend_ncol = 1
  )

  retangio_save_plot_pdf_png_safe(
    plot_obj = p + theme(legend.position = "none"),
    filename_base = file.path(outdir, paste0(stem, "_NOlegend")),
    width = retangio_width_main,
    height = retangio_height_main
  )

  p_leg <- retangio_make_manual_legend_plot(
    labels = color_levels,
    colors = palette_values,
    title_text = legend_title,
    ncol = 1
  )

  retangio_save_plot_pdf_png_safe(
    plot_obj = p_leg,
    filename_base = file.path(outdir, paste0(stem, "_legendOnly")),
    width = 5.5,
    height = max(4.5, length(color_levels) * 0.42)
  )

  if (make_split && split_col %in% colnames(obj@meta.data) && any(!is.na(obj@meta.data[[split_col]]))) {
    p_split <- DimPlot(
      object = obj,
      reduction = reduction,
      group.by = color_col,
      split.by = split_col,
      pt.size = retangio_pt_size,
      raster = FALSE,
      ncol = 3
    ) &
      retangio_theme_proj() &
      theme(
        legend.position = "none",
        strip.text = element_text(face = "bold")
      )

    retangio_save_plot_pdf_png_safe(
      plot_obj = p_split,
      filename_base = file.path(outdir, paste0(stem, "_split_by_condition_NOlegend")),
      width = 16,
      height = 8
    )
  }

  if (make_highlight && split_col %in% colnames(obj@meta.data) && any(!is.na(obj@meta.data[[split_col]]))) {
    present_levels <- split_levels[split_levels %in% unique(as.character(obj@meta.data[[split_col]]))]

    for (lev in present_levels) {
      p_hi <- retangio_plot_umap_highlight_subset_by_annotation(
        obj = obj,
        reduction = reduction,
        subset_col = split_col,
        subset_level = lev,
        color_col = color_col,
        color_levels = color_levels,
        palette_values = palette_values,
        xlim = lims$xlim,
        ylim = lims$ylim,
        title_text = paste0(title_text, " — highlight ", retangio_condition_labels[[lev]])
      ) + theme(legend.position = "none")

      retangio_save_plot_pdf_png_safe(
        plot_obj = p_hi,
        filename_base = file.path(outdir, paste0(stem, "_highlight_", lev, "_NOlegend")),
        width = retangio_width_main,
        height = retangio_height_main
      )
    }
  }
}
# ------------------------------------------------------------------
# Annotation abundance helpers
# ------------------------------------------------------------------

retangio_validate_annotation_levels <- function(values, allowed_levels, colname) {
  values <- unique(as.character(values))
  values <- values[!is.na(values)]

  missing_levels <- setdiff(values, allowed_levels)
  if (length(missing_levels) > 0) {
    stop(
      paste0(
        "Unexpected values found in ", colname, ": ",
        paste(missing_levels, collapse = ", ")
      )
    )
  }
}

retangio_compute_annotation_abundance <- function(
  obj,
  annotation_col,
  sample_col = "plot_condition",
  label_order
) {
  if (!annotation_col %in% colnames(obj@meta.data)) {
    stop(paste0(annotation_col, " not found in metadata."))
  }
  if (!sample_col %in% colnames(obj@meta.data)) {
    stop(paste0(sample_col, " not found in metadata."))
  }

  meta <- obj@meta.data

  annotation_vals <- as.character(meta[[annotation_col]])
  sample_vals <- as.character(meta[[sample_col]])

  retangio_validate_annotation_levels(
    values = annotation_vals,
    allowed_levels = label_order,
    colname = annotation_col
  )

  present_samples <- retangio_condition_levels[
    retangio_condition_levels %in% unique(sample_vals)
  ]

  if (length(present_samples) == 0) {
    stop("No valid plotting conditions found.")
  }

  keep_idx <- !is.na(annotation_vals) & !is.na(sample_vals) & sample_vals %in% present_samples
  meta_use <- data.frame(
    sample = factor(sample_vals[keep_idx], levels = present_samples),
    annotation = factor(annotation_vals[keep_idx], levels = label_order),
    stringsAsFactors = FALSE
  )

  # overall
  overall_tab <- table(
    factor(meta_use$annotation, levels = label_order)
  )

  overall_df <- data.frame(
    annotation = factor(label_order, levels = label_order),
    n_cells = as.integer(overall_tab[label_order]),
    stringsAsFactors = FALSE
  )
  overall_df$n_cells[is.na(overall_df$n_cells)] <- 0L
  overall_total <- sum(overall_df$n_cells)
  overall_df$fraction <- if (overall_total > 0) overall_df$n_cells / overall_total else 0
  overall_df$percent <- overall_df$fraction * 100

  # per sample
  per_sample_tab <- table(
    factor(meta_use$sample, levels = present_samples),
    factor(meta_use$annotation, levels = label_order)
  )

  per_sample_df <- as.data.frame(per_sample_tab, stringsAsFactors = FALSE)
  colnames(per_sample_df) <- c("sample", "annotation", "n_cells")

  per_sample_df$sample <- factor(as.character(per_sample_df$sample), levels = present_samples)
  per_sample_df$annotation <- factor(as.character(per_sample_df$annotation), levels = label_order)

  sample_totals <- per_sample_df %>%
    group_by(sample) %>%
    summarise(sample_total = sum(n_cells), .groups = "drop")

  per_sample_df <- per_sample_df %>%
    left_join(sample_totals, by = "sample") %>%
    mutate(
      fraction = ifelse(sample_total > 0, n_cells / sample_total, 0),
      percent = fraction * 100
    )

  list(
    overall = overall_df,
    per_sample = per_sample_df,
    present_samples = present_samples
  )
}

retangio_plot_annotation_abundance_overall <- function(
  overall_df,
  label_order,
  palette_values,
  title_text = "Annotated EC subtype abundance — overall"
) {
  overall_df <- overall_df %>%
    mutate(
      sample = factor("All samples", levels = "All samples"),
      annotation = factor(as.character(annotation), levels = label_order)
    )

  ggplot(overall_df, aes(x = sample, y = percent, fill = annotation)) +
    geom_col(width = 0.85) +
    scale_fill_manual(
      values = palette_values,
      breaks = label_order,
      limits = label_order,
      drop = FALSE,
      name = "Annotation"
    ) +
    scale_y_continuous(
      labels = scales::percent_format(scale = 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = title_text,
      x = NULL,
      y = "Relative abundance (%)"
    ) +
    retangio_theme_proj() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(face = "bold")
    )
}

retangio_plot_annotation_abundance_per_sample <- function(
  per_sample_df,
  sample_order,
  label_order,
  palette_values,
  title_text = "Annotated EC subtype abundance — by sample"
) {
  per_sample_df <- per_sample_df %>%
    mutate(
      sample = factor(as.character(sample), levels = sample_order),
      annotation = factor(as.character(annotation), levels = label_order)
    )

  ggplot(per_sample_df, aes(x = sample, y = percent, fill = annotation)) +
    geom_col(width = 0.85) +
    scale_x_discrete(labels = retangio_condition_labels[sample_order]) +
    scale_fill_manual(
      values = palette_values,
      breaks = label_order,
      limits = label_order,
      drop = FALSE,
      name = "Annotation"
    ) +
    scale_y_continuous(
      labels = scales::percent_format(scale = 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = title_text,
      x = NULL,
      y = "Relative abundance (%)"
    ) +
    retangio_theme_proj() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(face = "bold")
    )
}

retangio_plot_annotation_abundance_legend <- function(
  label_order,
  palette_values
) {
  legend_df <- data.frame(
    x = 1,
    y = seq_along(label_order),
    annotation = factor(label_order, levels = label_order)
  )

  p_leg <- ggplot(legend_df, aes(x = x, y = y, fill = annotation)) +
    geom_tile() +
    scale_fill_manual(
      values = palette_values,
      breaks = label_order,
      limits = label_order,
      drop = FALSE,
      name = "Annotation"
    ) +
    theme_void() +
    retangio_theme_proj() +
    theme(legend.position = "right")

  cowplot::plot_grid(retangio_extract_legend(p_leg), ncol = 1)
}

retangio_write_annotation_abundance_outputs <- function(
  obj,
  annotation_col,
  label_order,
  outdir_plots,
  outdir_results,
  plot_prefix,
  table_prefix,
  palette_values = NULL,
  title_prefix = "Annotated EC subtypes"
) {
  dir.create(outdir_plots, recursive = TRUE, showWarnings = FALSE)
  dir.create(outdir_results, recursive = TRUE, showWarnings = FALSE)

  if (is.null(palette_values)) {
    palette_values <- retangio_make_named_palette(label_order)
  }

  abundance_res <- retangio_compute_annotation_abundance(
    obj = obj,
    annotation_col = annotation_col,
    sample_col = "plot_condition",
    label_order = label_order
  )

  overall_df <- abundance_res$overall
  per_sample_df <- abundance_res$per_sample
  sample_order <- abundance_res$present_samples

  write.csv(
    overall_df,
    file = file.path(outdir_results, paste0(table_prefix, "_annotation_abundance_overall.csv")),
    row.names = FALSE
  )

  write.csv(
    per_sample_df,
    file = file.path(outdir_results, paste0(table_prefix, "_annotation_abundance_per_sample.csv")),
    row.names = FALSE
  )

  p_overall <- retangio_plot_annotation_abundance_overall(
    overall_df = overall_df,
    label_order = label_order,
    palette_values = palette_values,
    title_text = paste0(title_prefix, " — overall")
  )

  retangio_save_plot(
    plot_obj = p_overall,
    filename_base = file.path(outdir_plots, paste0(plot_prefix, "_annotation_abundance_overall_NOlegend")),
    width = 11,
    height = 8
  )

  p_per_sample <- retangio_plot_annotation_abundance_per_sample(
    per_sample_df = per_sample_df,
    sample_order = sample_order,
    label_order = label_order,
    palette_values = palette_values,
    title_text = paste0(title_prefix, " — by sample")
  )

  retangio_save_plot(
    plot_obj = p_per_sample,
    filename_base = file.path(outdir_plots, paste0(plot_prefix, "_annotation_abundance_per_sample_NOlegend")),
    width = 10,
    height = 8
  )

  p_legend <- retangio_plot_annotation_abundance_legend(
    label_order = label_order,
    palette_values = palette_values
  )

  retangio_save_plot(
    plot_obj = p_legend,
    filename_base = file.path(outdir_plots, paste0(plot_prefix, "_annotation_abundance_legendOnly")),
    width = 5,
    height = max(4, 0.45 * length(label_order))
  )
}
