#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

source("R/plotting_helpers.R")

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_rds <- "results/objects/10h_ec_after_soupx_harmony_annotated.rds"
  outdir <- "plots/10i_ec_after_soupx_harmony_annotated"

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 10i_plot_annotated_ec_after_soupx_harmony")

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  reduction_to_use <- if ("umap.harmony" %in% Reductions(seu)) "umap.harmony" else if ("umap" %in% Reductions(seu)) "umap" else stop("No Harmony UMAP reduction found.")
  seu <- retangio_add_plotting_metadata(seu)

  my_levels <- c(
    "Arterial x Capillary",
    "Capillary",
    "Venous x Capillary",
    "Tip cell",
    "Proliferative 1",
    "Proliferative 2",
    "Inflammatory",
    "EC 1",
    "EC 2",
    "Uncertain"
  )
  my_levels <- my_levels[my_levels %in% unique(as.character(seu$annotation_my))]
  my_pal <- retangio_make_named_palette(my_levels)

  proposed_levels <- c(
    "arterial_side_capillary_BRB_homeostatic",
    "capillary_BRB",
    "capillary_BRB_venous_biased",
    "tip_angiogenic",
    "proliferative_S_phase",
    "proliferative_G2M",
    "inflammatory_IFN_EC",
    "specialized_Aqp1_Six3_like_EC",
    "rare_specialized_EC",
    "atypical_capillary_BRB_stressed_mixed"
  )
  proposed_levels <- proposed_levels[proposed_levels %in% unique(as.character(seu$annotation_proposed))]
  proposed_pal <- retangio_make_named_palette(proposed_levels)

  retangio_plot_annotation_bundle(
    obj = seu,
    reduction = reduction_to_use,
    color_col = "annotation_my",
    color_levels = my_levels,
    palette_values = my_pal,
    outdir = outdir,
    stem = "10i_harmony_UMAP_annotation_my_labels",
    title_text = "UMAP — Harmony — annotated clusters (my labels)"
  )

  retangio_plot_annotation_bundle(
    obj = seu,
    reduction = reduction_to_use,
    color_col = "annotation_proposed",
    color_levels = proposed_levels,
    palette_values = proposed_pal,
    outdir = outdir,
    stem = "10i_harmony_UMAP_annotation_proposed_labels",
    title_text = "UMAP — Harmony — annotated clusters (proposed labels)"
  )

  message_ts("Stage 10i_plot_annotated_ec_after_soupx_harmony complete")
}

main()
