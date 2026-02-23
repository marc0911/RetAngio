source("R/00_utils.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

main <- function() {
  cfg <- read_config()
  outdir <- cfg$project$output_dir
  ensure_dir(outdir)
  set_seed_from_config(cfg)

  seu <- readRDS(file.path(outdir, "01_loaded_raw.rds"))
  cluster_col <- cfg$integration$contaminant_cluster_column
  contaminants <- as.character(cfg$integration$contaminants)
  sample_col <- cfg$input$sample_column
  dims <- parse_integer_sequence(cfg$integration$dims, "integration.dims")

  if (!cluster_col %in% colnames(seu@meta.data)) {
    fail_fast("Cluster column missing: {cluster_col}")
  }

  keep <- !(as.character(seu[[cluster_col]][, 1]) %in% contaminants)
  seu <- subset(seu, cells = colnames(seu)[keep])
  check_contaminants_removed(seu, cluster_col, contaminants)

  # Seurat v5 SCT+RPCA integration workflow:
  # split by sample -> SelectIntegrationFeatures -> RunPCA on selected SCT features -> IntegrateLayers.
  # Authoritative docs: Seurat v5 integration vignette
  # https://satijalab.org/seurat/articles/seurat5_integration.html
  obj_list <- SplitObject(seu, split.by = sample_col)
  features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = cfg$integration$nfeatures)

  seu <- merge(obj_list[[1]], y = obj_list[-1])
  VariableFeatures(seu) <- features

  seu <- RunPCA(seu, assay = cfg$integration$assay, features = features, npcs = max(dims), verbose = FALSE)

  seu[[cfg$integration$assay]] <- split(seu[[cfg$integration$assay]], f = seu[[sample_col]][, 1])
  seu <- IntegrateLayers(
    object = seu,
    method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.rpca",
    assay = cfg$integration$assay,
    dims = dims,
    verbose = FALSE
  )

  seu <- FindNeighbors(seu, reduction = "integrated.rpca", dims = dims)
  seu <- FindClusters(seu, resolution = cfg$integration$clustering_resolution)
  seu <- RunUMAP(seu, reduction = "integrated.rpca", dims = dims)

  # enforce final cluster column and Idents alignment
  new_cluster_col <- paste0("integrated_snn_res.", cfg$integration$clustering_resolution)
  if (!new_cluster_col %in% colnames(seu@meta.data)) {
    fail_fast("Expected cluster column missing after FindClusters: {new_cluster_col}")
  }
  seu[["final_cluster"]] <- seu[[new_cluster_col]]
  Idents(seu) <- "final_cluster"

  p <- DimPlot(seu, group.by = "final_cluster", label = TRUE) + ggtitle("Post-cleaning RPCA-integrated UMAP")
  ggsave(file.path(outdir, "02_umap_post_reintegration.png"), p, width = 8, height = 6, dpi = 300)

  saveRDS(seu, file.path(outdir, "02_clean_reintegrated.rds"))
  yaml::write_yaml(cfg, file.path(outdir, "config_used.yml"))
  log_info("02_remove_contaminants_and_reintegrate complete")
}

if (sys.nframe() == 0) {
  main()
}
