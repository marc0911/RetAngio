source("R/00_utils.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
})

# Method notes:
# - Seurat graph-based clustering/UMAP workflow follows Seurat v5 docs (authoritative Satija Lab docs).
# - ElbowPlot is used as a heuristic for PC selection (Seurat docs).

args <- parse_args()
cfg <- load_config(args$config)
set.seed(cfg$seed)
init_dirs(cfg)
script_name <- "02_remove_contaminants_and_recluster"

input_path <- file.path(cfg$output$root, cfg$output$snapshots_dir, "01_input_validated.rds")
if (!file.exists(input_path)) stop("Missing prior snapshot. Run 01 first.", call. = FALSE)
seu <- readRDS(input_path)

cluster_col <- cfg$contaminant_filter$cluster_column
clusters_to_remove <- as.character(cfg$contaminant_filter$clusters_to_remove)

pre_n <- ncol(seu)
keep_cells <- rownames(seu@meta.data)[!(as.character(seu@meta.data[[cluster_col]]) %in% clusters_to_remove)]
seu <- subset(seu, cells = keep_cells)
post_n <- ncol(seu)

if (post_n >= pre_n) stop("Contaminant filtering failed: cell count did not decrease.", call. = FALSE)

log_message(glue::glue("Removed {pre_n - post_n} cells from clusters: {paste(clusters_to_remove, collapse=',')}"), cfg, script_name)

DefaultAssay(seu) <- cfg$input_assays$default_assay
seu <- RunPCA(seu, npcs = cfg$recluster$npcs, verbose = FALSE, seed.use = cfg$seed)

elbow <- ElbowPlot(seu, ndims = cfg$recluster$npcs)
save_plot(elbow, "02_elbow_plot.png", cfg)

dims_grid <- cfg$recluster$dims_grid
res_grid <- cfg$recluster$resolution_grid
results <- vector("list", length(dims_grid) * length(res_grid))
idx <- 1

for (dims in dims_grid) {
  seu_dim <- FindNeighbors(seu, dims = 1:dims, k.param = cfg$recluster$k_param, verbose = FALSE)
  for (res in res_grid) {
    tmp <- FindClusters(seu_dim, resolution = res, verbose = FALSE)
    tmp <- RunUMAP(tmp, dims = 1:dims, seed.use = cfg$seed, verbose = FALSE)
    colname <- paste0("SCT_snn_res.", res)
    p <- DimPlot(tmp, group.by = colname, label = TRUE) +
      ggplot2::ggtitle(glue::glue("dims={dims}, res={res}"))
    file_tag <- str_replace_all(paste0("02_umap_dims", dims, "_res", res, ".png"), "\\.", "p")
    save_plot(p, file_tag, cfg)

    results[[idx]] <- tibble::tibble(
      dims = dims,
      resolution = res,
      n_clusters = dplyr::n_distinct(tmp@meta.data[[colname]]),
      mean_silhouette_proxy = NA_real_
    )
    idx <- idx + 1
  }
}

grid_summary <- dplyr::bind_rows(results)
summary_path <- file.path(cfg$output$root, cfg$output$tables_dir, "02_recluster_grid_summary.csv")
readr::write_csv(grid_summary, summary_path)

final_dims <- cfg$recluster$final_dims
final_res <- cfg$recluster$final_resolution
seu <- FindNeighbors(seu, dims = 1:final_dims, k.param = cfg$recluster$k_param, verbose = FALSE)
seu <- FindClusters(seu, resolution = final_res, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:final_dims, seed.use = cfg$seed, verbose = FALSE)

final_col <- paste0("SCT_snn_res.", final_res)
cfg$recluster$final_cluster_column <- final_col
Idents(seu) <- seu@meta.data[[final_col]]

final_plot <- DimPlot(seu, group.by = final_col, label = TRUE) + ggplot2::ggtitle("Final clustering")
save_plot(final_plot, "02_final_umap.png", cfg)

out <- file.path(cfg$output$root, cfg$output$snapshots_dir, "02_post_filter_reclustered.rds")
safe_save_rds(seu, out)

qc_tbl <- tibble::tibble(
  check = c("contaminants_removed", "idents_match_final_column"),
  pass = c(
    !any(as.character(seu@meta.data[[cluster_col]]) %in% clusters_to_remove),
    identical(as.character(Idents(seu)), as.character(seu@meta.data[[final_col]]))
  )
)
readr::write_csv(qc_tbl, file.path(cfg$output$root, cfg$output$reports_dir, "02_sanity_checks.csv"))

write_session_info(cfg, script_name)
log_message("Completed contaminant removal + reclustering", cfg, script_name)
