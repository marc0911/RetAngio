#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

source("R/plotting_helpers.R")

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_rds <- "results/objects/09h_ec_after_soupx_rpca_annotated.rds"

  output_rds <- "results/objects/09k_rpca_tipcell_resolution_sweep.rds"
  output_counts <- "results/review/09k_rpca_tipcell_counts_by_sample.csv"
  output_cluster_sizes <- "results/review/09k_rpca_tipcell_cluster_sizes_by_resolution.csv"
  plot_dir <- "plots/09k_rpca_tipcell_resolution_sweep"

  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/review", recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 09k_tipcell_subclustering_rpca")

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  if (!"annotation_my" %in% colnames(seu@meta.data)) {
    stop("annotation_my not found in metadata.")
  }

  tip_cells <- colnames(seu)[as.character(seu$annotation_my) == "Tip cell"]
  if (length(tip_cells) == 0) {
    stop("No cells found with annotation_my == 'Tip cell'.")
  }

  message_ts("Tip cells retained:", length(tip_cells))

  # ------------------------------------------------------------------
  # Build a fresh RNA-only object for robust local subclustering
  # ------------------------------------------------------------------
  DefaultAssay(seu) <- "RNA"
  message_ts("Joining RNA layers")
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])

  tip_meta <- seu@meta.data[tip_cells, , drop = FALSE]
  tip_counts <- GetAssayData(seu, assay = "RNA", layer = "counts")[, tip_cells, drop = FALSE]

  seu_tip <- CreateSeuratObject(
    counts = tip_counts,
    meta.data = tip_meta,
    project = "RPCA_tipcells"
  )

  seu_tip <- retangio_add_plotting_metadata(seu_tip)

  if ("sample_id" %in% colnames(seu_tip@meta.data)) {
    seu_tip$plot_condition_fixed <- factor(as.character(seu_tip$sample_id), levels = retangio_condition_levels)
  } else if ("orig.ident" %in% colnames(seu_tip@meta.data)) {
    seu_tip$plot_condition_fixed <- factor(as.character(seu_tip$orig.ident), levels = retangio_condition_levels)
  } else {
    stop("Neither sample_id nor orig.ident found in tip-cell object metadata.")
  }

  counts_by_sample <- seu_tip@meta.data %>%
    mutate(sample_plot = as.character(plot_condition_fixed)) %>%
    count(sample_plot, name = "n_cells") %>%
    mutate(sample_plot = factor(sample_plot, levels = retangio_condition_levels)) %>%
    arrange(sample_plot)

  write_csv(counts_by_sample, output_counts)
  message_ts("Saved tip-cell counts by sample:", output_counts)

  # ------------------------------------------------------------------
  # Fresh local clustering on RNA
  # ------------------------------------------------------------------
  DefaultAssay(seu_tip) <- "RNA"

  message_ts("Normalizing RNA")
  seu_tip <- NormalizeData(seu_tip, verbose = FALSE)

  message_ts("Finding variable features")
  seu_tip <- FindVariableFeatures(
    seu_tip,
    selection.method = "vst",
    nfeatures = 3000,
    verbose = FALSE
  )

  message_ts("Scaling data")
  seu_tip <- ScaleData(
    seu_tip,
    features = VariableFeatures(seu_tip),
    verbose = FALSE
  )

  message_ts("Running local PCA")
  seu_tip <- RunPCA(
    seu_tip,
    features = VariableFeatures(seu_tip),
    npcs = 30,
    reduction.name = "pca.tip",
    verbose = FALSE
  )

  message_ts("Running local neighbors")
  seu_tip <- FindNeighbors(
    seu_tip,
    reduction = "pca.tip",
    dims = 1:15,
    graph.name = c("tip_nn", "tip_snn"),
    verbose = FALSE
  )

  resolutions <- c(0.05, 0.10, 0.20, 0.40, 0.60)

  for (res in resolutions) {
    cluster_col <- paste0("tip_res_", gsub("\\.", "_", format(res, nsmall = 2)))
    message_ts("Clustering tip cells at resolution:", res, "->", cluster_col)

    seu_tip <- FindClusters(
      seu_tip,
      graph.name = "tip_snn",
      resolution = res,
      algorithm = 1,
      cluster.name = cluster_col,
      verbose = FALSE
    )
  }

  message_ts("Running local UMAP")
  seu_tip <- RunUMAP(
    seu_tip,
    reduction = "pca.tip",
    dims = 1:15,
    reduction.name = "umap.tip",
    reduction.key = "UMAPTIP_",
    verbose = FALSE
  )

  saveRDS(seu_tip, output_rds)
  message_ts("Saved tip-cell subcluster object:", output_rds)

  cluster_sizes_all <- bind_rows(lapply(resolutions, function(res) {
    cluster_col <- paste0("tip_res_", gsub("\\.", "_", format(res, nsmall = 2)))
    tibble(
      resolution = as.character(res),
      cluster = as.character(seu_tip@meta.data[[cluster_col]])
    ) %>%
      count(resolution, cluster, name = "n_cells", .drop = FALSE) %>%
      arrange(resolution, suppressWarnings(as.numeric(cluster)))
  }))

  write_csv(cluster_sizes_all, output_cluster_sizes)
  message_ts("Saved cluster sizes by resolution:", output_cluster_sizes)

  lims <- retangio_get_fixed_umap_limits(seu_tip, "umap.tip")

  # sample-colored UMAP
  df_cond <- retangio_make_umap_df(
    obj = seu_tip,
    reduction = "umap.tip",
    color_col = "plot_condition_fixed",
    levels_use = retangio_condition_levels,
    shuffle_seed = 1
  )

  p_cond <- retangio_plot_umap_discrete(
    df = df_cond,
    color_col = "plot_condition_fixed",
    palette_values = retangio_condition_palette,
    xlim = lims$xlim,
    ylim = lims$ylim,
    title_text = "RPCA tip-cell subset — by sample",
    legend_title = "Sample",
    legend_ncol = 1
  )

  retangio_save_plot(
    plot_obj = p_cond + theme(legend.position = "none"),
    filename_base = file.path(plot_dir, "09k_tipcell_plot_condition_NOlegend"),
    width = retangio_width_main,
    height = retangio_height_main
  )

  p_cond_leg <- retangio_make_manual_legend_plot(
    labels = retangio_condition_levels[retangio_condition_levels %in% unique(as.character(seu_tip$plot_condition_fixed))],
    colors = retangio_condition_palette,
    title_text = "Sample",
    ncol = 1
  )

  retangio_save_plot(
    plot_obj = p_cond_leg,
    filename_base = file.path(plot_dir, "09k_tipcell_plot_condition_legendOnly"),
    width = 4.5,
    height = 4.5
  )

  # subcluster UMAPs + highlight-by-condition
  for (res in resolutions) {
    cluster_col <- paste0("tip_res_", gsub("\\.", "_", format(res, nsmall = 2)))

    cluster_levels <- sort(unique(as.character(seu_tip@meta.data[[cluster_col]])))
    cluster_levels <- cluster_levels[!is.na(cluster_levels)]
    cluster_levels <- cluster_levels[order(suppressWarnings(as.numeric(cluster_levels)))]
    cluster_pal <- retangio_make_cluster_palette(cluster_levels)

    df_res <- retangio_make_umap_df(
      obj = seu_tip,
      reduction = "umap.tip",
      color_col = cluster_col,
      levels_use = cluster_levels,
      shuffle_seed = 1
    )

    p_res <- retangio_plot_umap_discrete(
      df = df_res,
      color_col = cluster_col,
      palette_values = cluster_pal,
      xlim = lims$xlim,
      ylim = lims$ylim,
      title_text = paste0("RPCA tip-cell subset — resolution ", format(res, nsmall = 2)),
      legend_title = "Subcluster",
      legend_ncol = 1
    )

    stem <- paste0("09k_tipcell_res_", gsub("\\.", "_", format(res, nsmall = 2)))

    retangio_save_plot(
      plot_obj = p_res + theme(legend.position = "none"),
      filename_base = file.path(plot_dir, paste0(stem, "_NOlegend")),
      width = retangio_width_main,
      height = retangio_height_main
    )

    p_leg <- retangio_make_manual_legend_plot(
      labels = cluster_levels,
      colors = cluster_pal,
      title_text = "Subcluster",
      ncol = 1
    )

    retangio_save_plot(
      plot_obj = p_leg,
      filename_base = file.path(plot_dir, paste0(stem, "_legendOnly")),
      width = 4.5,
      height = max(4.5, length(cluster_levels) * 0.45)
    )

    for (cond in retangio_condition_levels) {
      if (!(cond %in% unique(as.character(seu_tip$plot_condition_fixed)))) next

      df_hi <- df_res
      df_hi$plot_condition_fixed <- seu_tip@meta.data[rownames(df_hi), "plot_condition_fixed", drop = TRUE]
      df_hi$plot_condition_fixed <- factor(as.character(df_hi$plot_condition_fixed), levels = retangio_condition_levels)

      p_hi <- retangio_plot_umap_highlight(
        df = df_hi,
        group_col = "plot_condition_fixed",
        highlight_level = cond,
        palette_values = retangio_condition_palette,
        xlim = lims$xlim,
        ylim = lims$ylim,
        title_text = paste0(
          "RPCA tip-cell subset — resolution ",
          format(res, nsmall = 2),
          " — highlight ",
          retangio_condition_labels[[cond]]
        ),
        legend_title = "Sample"
      ) + theme(legend.position = "none")

      retangio_save_plot(
        plot_obj = p_hi,
        filename_base = file.path(plot_dir, paste0(stem, "_highlight_", cond, "_NOlegend")),
        width = retangio_width_main,
        height = retangio_height_main
      )
    }
  }

  message_ts("Stage 09k_tipcell_subclustering_rpca complete")
}

main()
