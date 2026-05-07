#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

repo_root <- "/work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo"
setwd(repo_root)

input_rds <- "results/objects/06_soupx_unintegrated_clustered_umap.rds"
input_markers_csv <- "results/markers/07a_unintegrated_after_soupx_findallmarkers.csv"

output_top20_csv <- "results/review/07b_cluster_review_top20_markers.csv"
output_panel_csv <- "results/review/07b_cluster_review_marker_panel_avgexp.csv"
output_sizes_csv <- "results/review/07b_cluster_review_cluster_sizes.csv"

dir.create("results/review", recursive = TRUE, showWarnings = FALSE)
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

log_message("Loading clustered object")
stop_if_missing(input_rds, "clustered object")
obj <- readRDS(input_rds)

log_message("Loading markers table")
stop_if_missing(input_markers_csv, "markers CSV")
markers <- read.csv(input_markers_csv, stringsAsFactors = FALSE, check.names = FALSE)

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

Idents(obj) <- "seurat_clusters"

log_message("Preparing top-20 marker review table")
top20 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(as.numeric(as.character(cluster)), desc(avg_log2FC), p_val_adj)

write.csv(top20, output_top20_csv, row.names = FALSE)

log_message("Preparing cluster size table")
cluster_sizes <- obj@meta.data %>%
  dplyr::count(seurat_clusters, name = "n_cells") %>%
  mutate(freq = n_cells / sum(n_cells)) %>%
  arrange(as.numeric(as.character(seurat_clusters)))

write.csv(cluster_sizes, output_sizes_csv, row.names = FALSE)

marker_panel <- c(
  "Cldn5", "Kdr", "Pecam1", "Klf2", "Emcn", "Cd34", "Ptprb", "Tek", "Vwf",
  "Apln", "Angpt2", "Adm", "Kit", "Esm1",
  "Gja5", "Gja4", "Efnb2", "Sox17", "Bmx", "Slco1a4",
  "Nr2f2", "Nrp2", "Vcam1", "Vwf",
  "Mki67", "Top2a", "Pclaf", "Birc5",
  "Pdgfrb", "Rgs5", "Cspg4", "Des", "Notch3",
  "Col1a1", "Col1a2", "Dcn", "Lum",
  "Rho", "Gnat1", "Pde6b", "Prph2", "Sag",
  "Snap25", "Syt1", "Rbfox3", "Elavl4",
  "Aqp4", "Rlbp1", "Slc1a3", "Glul",
  "Plp1", "Mbp", "Mog",
  "Lyz2", "Tyrobp", "C1qa", "Aif1",
  "Hbb-bs", "Hba-a1", "Hba-a2"
)

marker_panel <- unique(marker_panel)
marker_panel_present <- marker_panel[marker_panel %in% rownames(obj)]

log_message("Marker panel genes requested: ", length(marker_panel))
log_message("Marker panel genes present in object: ", length(marker_panel_present))

avg_exp <- AverageExpression(
  object = obj,
  assays = "RNA",
  features = marker_panel_present,
  group.by = "seurat_clusters",
  layer = "data",
  verbose = TRUE
)

panel_df <- as.data.frame(avg_exp$RNA)
panel_df$gene <- rownames(panel_df)
panel_df <- panel_df[, c("gene", setdiff(colnames(panel_df), "gene")), drop = FALSE]

write.csv(panel_df, output_panel_csv, row.names = FALSE)

log_message("Saved top-20 marker review table: ", output_top20_csv)
log_message("Saved marker panel average-expression table: ", output_panel_csv)
log_message("Saved cluster size table: ", output_sizes_csv)
log_message("Stage 07b completed successfully")
