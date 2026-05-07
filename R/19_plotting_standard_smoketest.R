#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

repo_root <- "/work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo"
setwd(repo_root)

source("R/plotting_helpers.R")

input_rds <- "results/objects/17_ec_integrated_clustered_umap.rds"
outdir <- "plots/19_plotting_standard_smoketest"

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

choose_group_col <- function(obj) {
  if ("plot_condition" %in% colnames(obj@meta.data) &&
      any(!is.na(obj@meta.data$plot_condition))) {
    return("plot_condition")
  }

  candidates <- c("condition", "group", "orig.ident", "sample_id")
  hit <- candidates[candidates %in% colnames(obj@meta.data)]
  if (length(hit) == 0) {
    stop("Could not find any grouping column among: plot_condition, condition, group, orig.ident, sample_id", call. = FALSE)
  }
  hit[[1]]
}

choose_umap_reduction <- function(obj) {
  red_names <- Reductions(obj)

  if ("umap" %in% red_names) {
    return("umap")
  }

  umap_hits <- grep("umap", red_names, ignore.case = TRUE, value = TRUE)

  if (length(umap_hits) == 1) {
    return(umap_hits[[1]])
  }

  if (length(umap_hits) > 1) {
    log_message("Multiple UMAP-like reductions found: ", paste(umap_hits, collapse = ", "))
    return(umap_hits[[1]])
  }

  stop(
    sprintf(
      "No UMAP reduction found in object. Available reductions: %s",
      paste(red_names, collapse = ", ")
    ),
    call. = FALSE
  )
}

log_message("Loading object: ", input_rds)
stop_if_missing(input_rds, "integrated clustered object")
obj <- readRDS(input_rds)

obj <- retangio_add_plotting_metadata(obj)

group_col <- choose_group_col(obj)
reduction_to_use <- choose_umap_reduction(obj)

if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
  obj$seurat_clusters <- as.character(Idents(obj))
}

if ("plot_condition_source" %in% colnames(obj@meta.data)) {
  log_message("plot_condition source: ", unique(as.character(obj$plot_condition_source))[1])
}

log_message("Using grouping column: ", group_col)
log_message("Using reduction: ", reduction_to_use)

df_group <- retangio_get_umap_df(
  obj = obj,
  reduction = reduction_to_use,
  color_col = group_col,
  order_schema = if (group_col == "plot_condition") "condition" else "auto",
  seed = 1234
)

lims <- retangio_get_xy_limits(df_group)

group_levels <- levels(df_group[[group_col]])

group_palette <- NULL
if (identical(group_levels, retangio_condition_levels)) {
  group_palette <- retangio_condition_palette
} else {
  group_palette <- retangio_make_named_palette(group_levels)
}

p_group <- retangio_plot_umap_discrete(
  df = df_group,
  color_col = group_col,
  palette_values = group_palette,
  title = "UMAP",
  subtitle = paste("Colored by", group_col),
  point_size = 0.10,
  point_alpha = 0.45,
  xlim = lims$xlim,
  ylim = lims$ylim,
  legend_title = group_col
)

retangio_save_plot_with_separate_legend(
  p = p_group,
  outdir = outdir,
  stem = paste0("UMAP_by_", group_col),
  width_main = 7,
  height_main = 6,
  width_legend = 4,
  height_legend = 3
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
  title = "UMAP",
  subtitle = "Colored by seurat_clusters",
  point_size = 0.10,
  point_alpha = 0.45,
  xlim = lims$xlim,
  ylim = lims$ylim,
  legend_title = "Cluster",
  base_size = 18
)

retangio_save_plot_with_separate_legend(
  p = p_cluster,
  outdir = outdir,
  stem = "UMAP_by_seurat_clusters",
  width_main = 7,
  height_main = 6,
  width_legend = 5,
  height_legend = 4
)

for (lev in group_levels) {
  p_hi <- retangio_plot_umap_highlight(
    df = df_group,
    group_col = group_col,
    highlight_level = lev,
    palette_values = group_palette,
    title = "UMAP",
    subtitle = paste("Highlight", lev),
    point_size = 0.10,
    highlight_alpha = 0.95,
    bg_size = 0.10,
    bg_alpha = 0.12,
    bg_color = "grey80",
    xlim = lims$xlim,
    ylim = lims$ylim,
    legend_title = group_col
  )

  safe_lev <- gsub("[^A-Za-z0-9_]+", "_", lev)

  ggsave(
    filename = file.path(outdir, paste0("UMAP_highlight_", safe_lev, "_NOlegend.pdf")),
    plot = p_hi + theme(legend.position = "none"),
    width = 7,
    height = 6,
    useDingbats = FALSE
  )
  ggsave(
    filename = file.path(outdir, paste0("UMAP_highlight_", safe_lev, "_NOlegend.png")),
    plot = p_hi + theme(legend.position = "none"),
    width = 7,
    height = 6,
    dpi = 300,
    bg = "white"
  )
}

log_message("Smoke-test plotting outputs written to: ", outdir)
log_message("Stage completed successfully")
