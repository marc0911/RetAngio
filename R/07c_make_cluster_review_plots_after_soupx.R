#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

repo_root <- "/work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo"
setwd(repo_root)

source("R/plotting_helpers.R")

input_rds <- "results/objects/06_soupx_unintegrated_clustered_umap.rds"
outdir <- "plots/07c_unintegrated_soupx_cluster_review"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

log_message <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", ts, paste(..., collapse = "")))
}

stop_if_missing <- function(path, label = NULL) {
  if (!file.exists(path) && !dir.exists(path)) {
    if (is.null(label)) {
      stop(sprintf("Required path does not exist: %s", path), call. = FALSE)
    } else {
      stop(sprintf("Required %s does not exist: %s", label, path), call. = FALSE)
    }
  }
}

choose_umap_reduction <- function(obj) {
  red_names <- Reductions(obj)

  if ("umap.unintegrated" %in% red_names) return("umap.unintegrated")
  if ("umap" %in% red_names) return("umap")

  umap_hits <- grep("umap", red_names, ignore.case = TRUE, value = TRUE)
  if (length(umap_hits) >= 1) return(umap_hits[[1]])

  stop(
    sprintf(
      "No UMAP reduction found in object. Available reductions: %s",
      paste(red_names, collapse = ", ")
    ),
    call. = FALSE
  )
}

save_feature_plot_single <- function(
  obj,
  feature,
  reduction_to_use,
  outdir,
  filename_base,
  min_cutoff,
  max_cutoff
) {
  p <- FeaturePlot(
    object = obj,
    features = feature,
    reduction = reduction_to_use,
    order = TRUE,
    min.cutoff = min_cutoff,
    max.cutoff = max_cutoff,
    pt.size = 0.1,
    combine = TRUE
  ) +
    retangio_theme(base_size = 16) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme(
      legend.position = "right"
    )

  retangio_save_plot(
    plot_obj = p,
    filename_base = file.path(outdir, filename_base),
    width = 6,
    height = 5
  )
}

save_violin_plot_single <- function(obj, feature, outdir, filename_base) {
  p <- VlnPlot(
    object = obj,
    features = feature,
    assay = "RNA",
    group.by = "seurat_clusters",
    pt.size = 0,
    combine = TRUE
  ) +
    retangio_theme(base_size = 16) +
    labs(x = "Cluster", y = "Expression level") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "none"
    )

  retangio_save_plot(
    plot_obj = p,
    filename_base = file.path(outdir, filename_base),
    width = 6,
    height = 5
  )
}


save_discrete_legend_plot <- function(labels, colors, title_text, outdir, filename_base, nrow_legend = NULL) {
  if (title_text == "Cluster") {
    labels <- as.character(sort(as.numeric(labels)))
  }

  df_leg <- data.frame(
    x = 1,
    y = seq_along(labels),
    lab = factor(labels, levels = labels)
  )

  if (is.null(nrow_legend)) {
    nrow_legend <- length(labels)
  }

  p <- ggplot(df_leg, aes(x = x, y = y, color = lab)) +
    geom_point(size = 5, alpha = 0) +
    scale_color_manual(
      values = colors[labels],
      breaks = labels,
      limits = labels,
      name = title_text
    ) +
    guides(color = guide_legend(
      ncol = 1,
      byrow = FALSE,
      override.aes = list(size = 6, alpha = 1)
    )) +
    retangio_theme(base_size = 16) +
    theme(
      legend.position = "right",
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )

  retangio_save_plot(
    plot_obj = p,
    filename_base = file.path(outdir, filename_base),
    width = 4.5,
    height = max(4, 0.42 * length(labels) + 1.5)
  )
}

get_family_q95_max <- function(obj, genes, assay = "RNA", layer = "data") {
  genes <- genes[genes %in% rownames(obj)]
  if (length(genes) == 0) return(1)

  mat <- GetAssayData(obj, assay = assay, layer = layer)[genes, , drop = FALSE]
  vals <- as.numeric(mat@x)

  if (length(vals) == 0) return(1)

  q <- as.numeric(stats::quantile(vals, probs = 0.95, na.rm = TRUE, names = FALSE))
  if (!is.finite(q) || q <= 0) q <- 1
  q
}

get_gene_q95_max <- function(obj, gene, assay = "RNA", layer = "data") {
  if (!gene %in% rownames(obj)) return(1)

  vec <- GetAssayData(obj, assay = assay, layer = layer)[gene, , drop = TRUE]
  vals <- as.numeric(vec)

  if (length(vals) == 0) return(1)

  q <- as.numeric(stats::quantile(vals, probs = 0.95, na.rm = TRUE, names = FALSE))
  if (!is.finite(q) || q <= 0) q <- 1
  q
}

log_message("Loading SoupX-corrected unintegrated clustered object")
stop_if_missing(input_rds, "SoupX-corrected unintegrated clustered object")
obj <- readRDS(input_rds)

if (!"RNA" %in% Assays(obj)) {
  stop("RNA assay not found in input object.", call. = FALSE)
}

DefaultAssay(obj) <- "RNA"

log_message("Joining RNA layers")
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

log_message("Normalizing joined RNA assay")
obj <- NormalizeData(
  object = obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
  obj$seurat_clusters <- as.character(Idents(obj))
}

obj <- retangio_add_plotting_metadata(obj)

reduction_to_use <- choose_umap_reduction(obj)
log_message("Using UMAP reduction: ", reduction_to_use)

group_col <- if ("plot_condition" %in% colnames(obj@meta.data) &&
                 any(!is.na(obj@meta.data$plot_condition))) {
  "plot_condition"
} else if ("sample_id" %in% colnames(obj@meta.data)) {
  "sample_id"
} else {
  "orig.ident"
}

log_message("Using sample/condition grouping column: ", group_col)

df_cond <- retangio_get_umap_df(
  obj = obj,
  reduction = reduction_to_use,
  color_col = group_col,
  order_schema = if (group_col == "plot_condition") "condition" else "alphabetical",
  seed = 1234
)

lims <- retangio_get_xy_limits(df_cond)

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
  filename_base = file.path(outdir, "07c_unintegrated_soupx_UMAP_by_condition_NOlegend"),
  width = 7,
  height = 6
)

save_discrete_legend_plot(
  labels = retangio_condition_levels,
  colors = retangio_condition_palette,
  title_text = "Condition",
  outdir = outdir,
  filename_base = "07c_unintegrated_soupx_UMAP_by_condition_legendOnly"
)

df_cluster <- retangio_get_umap_df(
  obj = obj,
  reduction = reduction_to_use,
  color_col = "seurat_clusters",
  order_schema = "alphabetical",
  seed = 1234
)

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
  filename_base = file.path(outdir, "07c_unintegrated_soupx_UMAP_by_cluster_NOlegend"),
  width = 7,
  height = 6
)

save_discrete_legend_plot(
  labels = as.character(sort(as.numeric(levels(df_cluster$seurat_clusters)))),
  colors = cluster_palette,
  title_text = "Cluster",
  outdir = outdir,
  filename_base = "07c_unintegrated_soupx_UMAP_by_cluster_legendOnly"
)

marker_families <- list(
  EC_core = c("Kdr", "Ptprb", "Pecam1", "Cldn5", "Emcn", "Cd34", "Tek", "Vwf"),
  Tip_activated = c("Apln", "Angpt2", "Esm1", "Adm", "Pgf"),
  Arterial_specialized = c("Sox17", "Gja4", "Bmx", "Efnb2", "Slco1a4"),
  Proliferation = c("Mki67", "Top2a", "Pclaf", "Birc5"),
  Mural = c("Rgs5", "Cspg4", "Pdgfrb", "Des", "Notch3"),
  Neuronal_photoreceptor = c("Snap25", "Syt1", "Rho", "Pde6b", "Prph2"),
  Glial = c("Rlbp1", "Slc1a3", "Glul", "Aqp4")
)

marker_families <- lapply(marker_families, function(x) x[x %in% rownames(obj)])
marker_families <- marker_families[vapply(marker_families, length, integer(1)) > 0]

family_q95_max <- lapply(marker_families, function(genes) {
  get_family_q95_max(obj, genes, assay = "RNA", layer = "data")
})

for (fam in names(marker_families)) {
  log_message("Processing marker family: ", fam)
  fam_genes <- marker_families[[fam]]
  fam_max <- family_q95_max[[fam]]

  log_message("Family ", fam, " comparable max (q95 pooled) = ", signif(fam_max, 4))

  for (gene in fam_genes) {
    gene_max <- get_gene_q95_max(obj, gene, assay = "RNA", layer = "data")

    log_message(
      "Gene ", gene,
      " visibility max (gene q95) = ", signif(gene_max, 4),
      " | comparable max (family q95) = ", signif(fam_max, 4)
    )

    save_feature_plot_single(
      obj = obj,
      feature = gene,
      reduction_to_use = reduction_to_use,
      outdir = outdir,
      filename_base = paste0("07c_unintegrated_soupx_FeaturePlot_visibility_", fam, "_", gene),
      min_cutoff = 0,
      max_cutoff = gene_max
    )

    save_feature_plot_single(
      obj = obj,
      feature = gene,
      reduction_to_use = reduction_to_use,
      outdir = outdir,
      filename_base = paste0("07c_unintegrated_soupx_FeaturePlot_comparable_", fam, "_", gene),
      min_cutoff = 0,
      max_cutoff = fam_max
    )

    save_violin_plot_single(
      obj = obj,
      feature = gene,
      outdir = outdir,
      filename_base = paste0("07c_unintegrated_soupx_VlnPlot_", fam, "_", gene)
    )
  }
}

dotplot_features <- marker_families

log_message("Generating dot plot")
p_dot <- DotPlot(
  object = obj,
  features = dotplot_features,
  assay = "RNA",
  group.by = "seurat_clusters"
) +
  retangio_theme(base_size = 12) +
  labs(x = "Markers", y = "Cluster") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

retangio_save_plot(
  plot_obj = p_dot,
  filename_base = file.path(outdir, "07c_unintegrated_soupx_DotPlot_cluster_identification"),
  width = 16,
  height = 8
)

log_message("Stage 07c completed successfully")
