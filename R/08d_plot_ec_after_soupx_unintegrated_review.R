#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

source("R/plotting_helpers.R")

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_rds <- "results/objects/08c_ec_after_soupx_unintegrated_clustered_umap.rds"
  outdir <- "plots/08d_ec_after_soupx_unintegrated_review"

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 08d_plot_ec_after_soupx_unintegrated_review")

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  seu <- retangio_add_plotting_metadata(seu)

  reduction_to_use <- if ("umap.unintegrated" %in% Reductions(seu)) {
    "umap.unintegrated"
  } else if ("umap" %in% Reductions(seu)) {
    "umap"
  } else {
    stop("No UMAP reduction found in object.")
  }

  message_ts("Using reduction:", reduction_to_use)

  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    seu$seurat_clusters <- as.character(Idents(seu))
  }

  df_cluster <- retangio_get_umap_df(
    obj = seu,
    reduction = reduction_to_use,
    color_col = "seurat_clusters",
    order_schema = "alphabetical",
    seed = 1234
  )

  lims <- retangio_get_xy_limits(df_cluster)
  cluster_palette <- retangio_make_named_palette(levels(df_cluster$seurat_clusters))

  p_cluster <- retangio_plot_umap_discrete(
    df = df_cluster,
    color_col = "seurat_clusters",
    palette_values = cluster_palette,
    point_size = 0.10,
    point_alpha = 0.45,
    xlim = lims$xlim,
    ylim = lims$ylim,
    legend_title = "Cluster",
    base_size = 18
  )

  retangio_save_plot(
    plot_obj = p_cluster + theme(legend.position = "none"),
    filename_base = file.path(outdir, "08d_ec_after_soupx_unintegrated_UMAP_by_cluster_NOlegend"),
    width = 7,
    height = 6
  )

  group_col <- if ("plot_condition" %in% colnames(seu@meta.data) &&
                   any(!is.na(seu@meta.data$plot_condition))) {
    "plot_condition"
  } else if ("sample_id" %in% colnames(seu@meta.data)) {
    "sample_id"
  } else {
    "orig.ident"
  }

  df_cond <- retangio_get_umap_df(
    obj = seu,
    reduction = reduction_to_use,
    color_col = group_col,
    order_schema = if (group_col == "plot_condition") "condition" else "alphabetical",
    seed = 1234
  )

  group_palette <- if (group_col == "plot_condition") {
    retangio_condition_palette
  } else {
    retangio_make_named_palette(levels(df_cond[[group_col]]))
  }

  p_cond <- retangio_plot_umap_discrete(
    df = df_cond,
    color_col = group_col,
    palette_values = group_palette,
    point_size = 0.10,
    point_alpha = 0.45,
    xlim = lims$xlim,
    ylim = lims$ylim,
    legend_title = if (group_col == "plot_condition") "Condition" else group_col,
    base_size = 18
  )

  retangio_save_plot(
    plot_obj = p_cond + theme(legend.position = "none"),
    filename_base = file.path(outdir, "08d_ec_after_soupx_unintegrated_UMAP_by_condition_NOlegend"),
    width = 7,
    height = 6
  )

  message_ts("Stage 08d_plot_ec_after_soupx_unintegrated_review complete")
}

main()
