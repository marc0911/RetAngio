#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

source("scripts/00_retangio_plotting_helpers.R")

main <- function() {
  input_rds <- "results/objects/09o_rpca_tipcell_ctrl_only_resolution_sweep.rds"
  outdir <- "plots/10x_presentation_tipcells_ctrl_only"
  out_csv <- "results/abundance/10x_tipcells_ctrl_only_cluster_relative_abundance_per_sample.csv"
  cluster_col <- "tip_ctrl_res_0_20"

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create("results/abundance", recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 10x_presentation_plots_tipcells_ctrl_only")

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  seu <- retangio_add_plotting_metadata(seu)

  reduction_to_use <- retangio_detect_reduction(
    seu,
    preferred = c("umap.tip.ctrl", "umap.tip", "umap")
  )

  lims <- retangio_get_fixed_umap_limits(seu, reduction_to_use)

  sample_levels_present <- retangio_condition_levels[
    retangio_condition_levels %in% unique(as.character(seu$plot_condition))
  ]

  df_sample <- retangio_make_umap_df(
    obj = seu,
    reduction = reduction_to_use,
    color_col = "plot_condition",
    levels_use = sample_levels_present,
    shuffle_seed = 1
  )

  p_sample <- retangio_plot_umap_discrete(
    df = df_sample,
    color_col = "plot_condition",
    palette_values = retangio_condition_palette[sample_levels_present],
    xlim = lims$xlim,
    ylim = lims$ylim,
    title_text = "RPCA control-only tip-cell subset — colored by sample",
    legend_title = "Sample"
  )

  retangio_save_plot(
    plot_obj = p_sample + theme(legend.position = "none"),
    filename_base = file.path(outdir, "10x_tipcells_ctrl_only_umap_by_sample_NOlegend"),
    width = 10,
    height = 8
  )

  retangio_save_plot(
    plot_obj = retangio_extract_legend_plot(p_sample),
    filename_base = file.path(outdir, "10x_tipcells_ctrl_only_umap_by_sample_legendOnly"),
    width = 4,
    height = 5
  )

  abundance_df <- retangio_make_abundance_df(seu, cluster_col = cluster_col)
  write.csv(abundance_df, out_csv, row.names = FALSE)

  p_abund <- retangio_plot_abundance(
    df = abundance_df,
    title_text = "Control-only tip-cell local clusters — relative abundance by sample",
    legend_title = "Tip cluster"
  )

  retangio_save_plot(
    plot_obj = p_abund + theme(legend.position = "none"),
    filename_base = file.path(outdir, "10x_tipcells_ctrl_only_cluster_relative_abundance_per_sample_NOlegend"),
    width = 10,
    height = 8
  )

  retangio_save_plot(
    plot_obj = retangio_extract_legend_plot(p_abund),
    filename_base = file.path(outdir, "10x_tipcells_ctrl_only_cluster_relative_abundance_per_sample_legendOnly"),
    width = 4,
    height = 6
  )

  message_ts("Saved abundance table:", out_csv)
  message_ts("Stage 10x_presentation_plots_tipcells_ctrl_only complete")
}

main()
