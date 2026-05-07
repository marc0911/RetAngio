source("R/00_utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tibble)
  library(future)
})

main <- function() {
  log_message("Starting stage 09_rpca_fasttrack_v5")

  cfg <- read_config()
  outdir <- cfg$project$output_dir
  sample_col <- cfg$analysis$sample_column

  ensure_dir(file.path(outdir, "objects"))
  ensure_dir(file.path(outdir, "integration"))
  ensure_dir(file.path("plots", "integration"))

  in_rds <- file.path(outdir, "objects", "04_qc_doublet_filtered.rds")
  if (!file.exists(in_rds)) {
    stop("Missing input object: ", in_rds)
  }

  seu <- readRDS(in_rds)
  log_message("Loaded object: ", in_rds)

  if (!sample_col %in% colnames(seu@meta.data)) {
    stop("Sample column not found in metadata: ", sample_col)
  }

  # Memory-safe behavior on cluster
  future::plan("sequential")
  options(future.globals.maxSize = 8 * 1024^3)

  DefaultAssay(seu) <- "RNA"

  if ("RNA" %in% names(seu@assays)) {
    log_message("Joining RNA layers before split")
    seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  }

  # 1) Split by sample
  seu_list <- SplitObject(seu, split.by = sample_col)
  log_message("Split object into ", length(seu_list), " sample-level objects")

  # 2) SCTransform per sample
  for (sid in names(seu_list)) {
    log_message("Running SCTransform for sample: ", sid)
    DefaultAssay(seu_list[[sid]]) <- "RNA"
    seu_list[[sid]] <- SCTransform(
      seu_list[[sid]],
      assay = "RNA",
      new.assay.name = "SCT",
      vst.flavor = "v2",
      verbose = FALSE
    )
  }

  # 3) Shared feature set
  nfeatures_use <- 2000
  features <- SelectIntegrationFeatures(
    object.list = seu_list,
    nfeatures = nfeatures_use
  )
  log_message("Selected integration features: ", length(features))

  # 4) Merge SCT objects (still NOT integrated)
  seu_pre <- merge(
    x = seu_list[[1]],
    y = seu_list[-1]
  )
  DefaultAssay(seu_pre) <- "SCT"

  # 5) Explicit variable features + PCA on merged pre-integration object
  VariableFeatures(seu_pre) <- features
  seu_pre <- RunPCA(
    seu_pre,
    assay = "SCT",
    features = features,
    npcs = 30,
    verbose = FALSE
  )
  log_message("Computed PCA on merged pre-integration object")

  # Quick pre-integration UMAP for comparison
  dims_use <- 1:12
  seu_pre <- FindNeighbors(
    seu_pre,
    reduction = "pca",
    dims = dims_use,
    verbose = FALSE
  )
  seu_pre <- FindClusters(
    seu_pre,
    resolution = 0.4,
    verbose = FALSE
  )
  seu_pre <- RunUMAP(
    seu_pre,
    reduction = "pca",
    dims = dims_use,
    reduction.name = "umap.pre",
    reduction.key = "UMAPPRE_",
    verbose = FALSE
  )

  saveRDS(seu_pre, file.path(outdir, "objects", "09_preintegration_sct_merged.rds"))
  log_message("Saved pre-integration object: results/objects/09_preintegration_sct_merged.rds")

  # 6) Split SCT assay by sample for Seurat v5 integration
  seu_pre[["SCT"]] <- split(seu_pre[["SCT"]], f = seu_pre[[sample_col]][,1])
  log_message("Split SCT assay into layers by sample for IntegrateLayers")

  # 7) RPCA integration
  seu_int <- IntegrateLayers(
    object = seu_pre,
    method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.rpca",
    assay = "SCT",
    verbose = TRUE
  )

  # 8) Clustering / UMAP on integrated reduction
  DefaultAssay(seu_int) <- "SCT"

  seu_int <- FindNeighbors(
    seu_int,
    reduction = "integrated.rpca",
    dims = dims_use,
    graph.name = "integrated_snn",
    verbose = FALSE
  )

  seu_int <- FindClusters(
    seu_int,
    graph.name = "integrated_snn",
    resolution = 0.4,
    algorithm = 1,
    verbose = FALSE
  )

  seu_int <- RunUMAP(
    seu_int,
    reduction = "integrated.rpca",
    dims = dims_use,
    reduction.name = "umap.integrated",
    reduction.key = "UMAPI_",
    verbose = FALSE
  )

  saveRDS(seu_int, file.path(outdir, "objects", "09_rpca_integrated_fasttrack.rds"))
  log_message("Saved integrated object: results/objects/09_rpca_integrated_fasttrack.rds")

  write_csv(
    tibble(feature = features),
    file.path(outdir, "integration", "09_rpca_features_used.csv")
  )

  cluster_col <- "integrated_snn_res.0.4"
  if (!cluster_col %in% colnames(seu_int@meta.data)) {
    stop("Missing integrated cluster column: ", cluster_col)
  }

  cluster_tbl <- seu_int@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    transmute(
      cell_id = cell_id,
      sample_id = .data[[sample_col]],
      cluster = as.character(.data[[cluster_col]])
    ) %>%
    count(cluster, sample_id, name = "n") %>%
    arrange(as.numeric(cluster), sample_id)

  write_csv(cluster_tbl, file.path(outdir, "integration", "09_rpca_cluster_counts_by_sample.csv"))

  cluster_totals <- seu_int@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    transmute(cluster = as.character(.data[[cluster_col]])) %>%
    count(cluster, name = "n") %>%
    arrange(as.numeric(cluster))

  write_csv(cluster_totals, file.path(outdir, "integration", "09_rpca_cluster_counts_total.csv"))

  # Plots
  p_pre_sample <- DimPlot(
    seu_pre,
    reduction = "umap.pre",
    group.by = sample_col
  ) + ggtitle("Pre-integration UMAP by sample")

  ggsave(
    filename = file.path("plots", "integration", "09_preintegration_umap_by_sample.png"),
    plot = p_pre_sample,
    width = 8,
    height = 6,
    dpi = 300
  )

  p_pre_cluster <- DimPlot(
    seu_pre,
    reduction = "umap.pre",
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE
  ) + ggtitle("Pre-integration UMAP clusters")

  ggsave(
    filename = file.path("plots", "integration", "09_preintegration_umap_clusters.png"),
    plot = p_pre_cluster,
    width = 8,
    height = 6,
    dpi = 300
  )

  p_int_sample <- DimPlot(
    seu_int,
    reduction = "umap.integrated",
    group.by = sample_col
  ) + ggtitle("RPCA integrated UMAP by sample")

  ggsave(
    filename = file.path("plots", "integration", "09_rpca_umap_by_sample.png"),
    plot = p_int_sample,
    width = 8,
    height = 6,
    dpi = 300
  )

  p_int_cluster <- DimPlot(
    seu_int,
    reduction = "umap.integrated",
    group.by = cluster_col,
    label = TRUE,
    repel = TRUE
  ) + ggtitle("RPCA integrated UMAP (res 0.4)")

  ggsave(
    filename = file.path("plots", "integration", "09_rpca_umap_clusters_res0.4.png"),
    plot = p_int_cluster,
    width = 8,
    height = 6,
    dpi = 300
  )

  p_int_split <- DimPlot(
    seu_int,
    reduction = "umap.integrated",
    group.by = cluster_col,
    split.by = sample_col,
    ncol = 3
  ) + ggtitle("RPCA integrated UMAP split by sample")

  ggsave(
    filename = file.path("plots", "integration", "09_rpca_umap_clusters_split_samples.png"),
    plot = p_int_split,
    width = 14,
    height = 8,
    dpi = 300
  )

  log_message("Stage 09_rpca_fasttrack_v5 complete")
}

main()
