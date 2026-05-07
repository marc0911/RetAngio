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

save_plot_safe <- function(plot_obj, filename_base, width, height, dpi = 300) {
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

main <- function() {
  input_rds <- "results/objects/09o_rpca_tipcell_ctrl_only_resolution_sweep.rds"

  out_rds <- "results/objects/09p_rpca_tipcell_ctrl_only_zarkada_scored.rds"
  out_sig <- "results/signatures/09p_rpca_tipcell_ctrl_only_zarkada_signatures_used.csv"
  out_cell <- "results/signatures/09p_rpca_tipcell_ctrl_only_zarkada_cell_scores.csv"
  out_cluster_mean <- "results/signatures/09p_rpca_tipcell_ctrl_only_zarkada_cluster_mean_scores.csv"
  out_cluster_sizes <- "results/review/09p_rpca_tipcell_ctrl_only_cluster_sizes_res_0_20.csv"
  plot_dir <- "plots/09p_rpca_tipcell_ctrl_only_zarkada"

  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/signatures", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/review", recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 09p_score_zarkada_tip_programs_rpca_ctrl_only")

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  cluster_col <- "tip_ctrl_res_0_20"
  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop(paste0(cluster_col, " not found in metadata."))
  }
  if (!"umap.tip.ctrl" %in% Reductions(seu)) {
    stop("umap.tip.ctrl reduction not found.")
  }

  seu <- retangio_add_plotting_metadata(seu)

  s_tip_genes_raw <- c(
    "Esm1", "Apln", "Angpt2", "Adm", "Kcne3",
    "Trp53i11", "Prnd", "Pgf", "Thbs1", "Edn1"
  )

  d_tip_genes_raw <- c(
    "Cldn5", "Cd34", "Cxcr4", "Dll4", "Pdgfb",
    "Efnb2", "Igf2", "Igf1r", "Kdr", "Flt4"
  )

  row_genes <- rownames(seu)

  map_present <- function(genes, universe) {
    unique(genes[genes %in% universe])
  }

  s_tip_genes <- map_present(s_tip_genes_raw, row_genes)
  d_tip_genes <- map_present(d_tip_genes_raw, row_genes)

  if (length(s_tip_genes) < 3) stop("Too few S-tip genes found in object.")
  if (length(d_tip_genes) < 3) stop("Too few D-tip genes found in object.")

  sig_df <- bind_rows(
    data.frame(signature = "zarkada_S_tip", gene = s_tip_genes, stringsAsFactors = FALSE),
    data.frame(signature = "zarkada_D_tip", gene = d_tip_genes, stringsAsFactors = FALSE)
  )
  write_csv(sig_df, out_sig)
  message_ts("Saved signature genes used:", out_sig)

  DefaultAssay(seu) <- "RNA"
  seu <- NormalizeData(seu, verbose = FALSE)

  message_ts("Scoring Zarkada S-tip program")
  seu <- AddModuleScore(
    object = seu,
    features = list(s_tip_genes),
    name = "zarkada_S_tip",
    assay = "RNA",
    search = FALSE
  )

  message_ts("Scoring Zarkada D-tip program")
  seu <- AddModuleScore(
    object = seu,
    features = list(d_tip_genes),
    name = "zarkada_D_tip",
    assay = "RNA",
    search = FALSE
  )

  s_col <- "zarkada_S_tip1"
  d_col <- "zarkada_D_tip1"

  seu$tip_ctrl_cluster_0_20 <- factor(
    as.character(seu@meta.data[[cluster_col]]),
    levels = sort(unique(as.character(seu@meta.data[[cluster_col]])), na.last = TRUE)
  )

  sample_col <- if ("sample_id" %in% colnames(seu@meta.data)) "sample_id" else "orig.ident"

  cell_df <- seu@meta.data %>%
    tibble::rownames_to_column("cell_barcode") %>%
    transmute(
      cell_barcode = cell_barcode,
      sample_id = as.character(.data[[sample_col]]),
      plot_condition = as.character(plot_condition_fixed),
      tip_ctrl_cluster_0_20 = as.character(tip_ctrl_cluster_0_20),
      zarkada_S_tip_score = .data[[s_col]],
      zarkada_D_tip_score = .data[[d_col]]
    )

  write_csv(cell_df, out_cell)
  message_ts("Saved cell-level scores:", out_cell)

  cluster_mean_df <- cell_df %>%
    group_by(tip_ctrl_cluster_0_20) %>%
    summarise(
      n_cells = n(),
      mean_S_tip_score = mean(zarkada_S_tip_score, na.rm = TRUE),
      mean_D_tip_score = mean(zarkada_D_tip_score, na.rm = TRUE),
      mean_D_minus_S = mean(zarkada_D_tip_score - zarkada_S_tip_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(suppressWarnings(as.numeric(tip_ctrl_cluster_0_20)))

  write_csv(cluster_mean_df, out_cluster_mean)
  message_ts("Saved cluster mean scores:", out_cluster_mean)

  cluster_sizes_df <- cell_df %>%
    count(tip_ctrl_cluster_0_20, name = "n_cells") %>%
    arrange(suppressWarnings(as.numeric(tip_ctrl_cluster_0_20)))

  write_csv(cluster_sizes_df, out_cluster_sizes)
  message_ts("Saved cluster sizes:", out_cluster_sizes)

  saveRDS(seu, out_rds)
  message_ts("Saved scored object:", out_rds)

  lims <- retangio_get_fixed_umap_limits(seu, "umap.tip.ctrl")

  cluster_levels <- levels(seu$tip_ctrl_cluster_0_20)
  cluster_pal <- retangio_make_cluster_palette(cluster_levels)

  df_cluster <- retangio_make_umap_df(
    obj = seu,
    reduction = "umap.tip.ctrl",
    color_col = "tip_ctrl_cluster_0_20",
    levels_use = cluster_levels,
    shuffle_seed = 1
  )

  p_cluster <- retangio_plot_umap_discrete(
    df = df_cluster,
    color_col = "tip_ctrl_cluster_0_20",
    palette_values = cluster_pal,
    xlim = lims$xlim,
    ylim = lims$ylim,
    title_text = "RPCA control tip-cell subset — local clusters (res 0.20)",
    legend_title = "Tip cluster",
    legend_ncol = 1
  )

  save_plot_safe(
    p_cluster + theme(legend.position = "none"),
    file.path(plot_dir, "09p_tip_ctrl_cluster_res_0_20_NOlegend"),
    width = retangio_width_main,
    height = retangio_height_main
  )

  p_cluster_leg <- retangio_make_manual_legend_plot(
    labels = cluster_levels,
    colors = cluster_pal,
    title_text = "Tip cluster",
    ncol = 1
  )

  save_plot_safe(
    p_cluster_leg,
    file.path(plot_dir, "09p_tip_ctrl_cluster_res_0_20_legendOnly"),
    width = 4.5,
    height = max(4.5, length(cluster_levels) * 0.45)
  )

  p_s <- FeaturePlot(
    object = seu,
    reduction = "umap.tip.ctrl",
    features = s_col,
    pt.size = retangio_pt_size,
    raster = FALSE
  ) +
    ggtitle("Control-only Zarkada S-tip score") +
    retangio_theme_proj()

  save_plot_safe(
    p_s,
    file.path(plot_dir, "09p_zarkada_S_tip_score_featureplot"),
    width = retangio_width_main,
    height = retangio_height_main
  )

  p_d <- FeaturePlot(
    object = seu,
    reduction = "umap.tip.ctrl",
    features = d_col,
    pt.size = retangio_pt_size,
    raster = FALSE
  ) +
    ggtitle("Control-only Zarkada D-tip score") +
    retangio_theme_proj()

  save_plot_safe(
    p_d,
    file.path(plot_dir, "09p_zarkada_D_tip_score_featureplot"),
    width = retangio_width_main,
    height = retangio_height_main
  )

  p_s_vln <- VlnPlot(
    object = seu,
    features = s_col,
    group.by = "tip_ctrl_cluster_0_20",
    pt.size = 0,
    assay = "RNA"
  ) +
    ggtitle("Control-only Zarkada S-tip score by local tip cluster") +
    xlab("tip_ctrl_res_0_20") +
    ylab("Module score") +
    retangio_theme_proj()

  save_plot_safe(
    p_s_vln,
    file.path(plot_dir, "09p_zarkada_S_tip_score_violin_by_cluster"),
    width = 9,
    height = 6
  )

  p_d_vln <- VlnPlot(
    object = seu,
    features = d_col,
    group.by = "tip_ctrl_cluster_0_20",
    pt.size = 0,
    assay = "RNA"
  ) +
    ggtitle("Control-only Zarkada D-tip score by local tip cluster") +
    xlab("tip_ctrl_res_0_20") +
    ylab("Module score") +
    retangio_theme_proj()

  save_plot_safe(
    p_d_vln,
    file.path(plot_dir, "09p_zarkada_D_tip_score_violin_by_cluster"),
    width = 9,
    height = 6
  )

  dot_features <- unique(c(s_tip_genes, d_tip_genes))

  p_dot <- DotPlot(
    object = seu,
    features = dot_features,
    group.by = "tip_ctrl_cluster_0_20",
    assay = "RNA"
  ) +
    RotatedAxis() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    ggtitle("Control-only Zarkada S-tip / D-tip genes by local tip cluster") +
    retangio_theme_proj()

  save_plot_safe(
    p_dot,
    file.path(plot_dir, "09p_zarkada_tip_gene_dotplot_by_cluster"),
    width = 12,
    height = 7
  )

  message_ts("Stage 09p_score_zarkada_tip_programs_rpca_ctrl_only complete")
}

main()
