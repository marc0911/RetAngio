source("R/00_utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tibble)
  library(patchwork)
  library(future)
})

main <- function() {
  log_message("Starting stage 10_review_contaminants_unintegrated")

  cfg <- read_config()
  outdir <- cfg$project$output_dir
  sample_col <- cfg$analysis$sample_column

  ensure_dir(file.path(outdir, "markers"))
  ensure_dir(file.path("plots", "markers"))

  in_rds <- file.path(outdir, "objects", "08_unintegrated_clustered_umap.rds")
  if (!file.exists(in_rds)) stop("Missing input object: ", in_rds)

  seu <- readRDS(in_rds)
  log_message("Loaded object: ", in_rds)

  cluster_col <- "unintegrated_snn_res.0.4"
  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop("Missing cluster column: ", cluster_col)
  }

  # Memory-safe DE behavior
  future::plan("sequential")
  options(future.globals.maxSize = 8 * 1024^3)

  # Use RNA for marker detection / annotation
  if (!"RNA" %in% names(seu@assays)) {
    stop("RNA assay not found in object")
  }

  DefaultAssay(seu) <- "RNA"

  # Seurat v5 safeguard: join layers, then create RNA data layer
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)

  Idents(seu) <- seu@meta.data[[cluster_col]]

  # Report cluster sizes first
  cluster_sizes <- table(Idents(seu))
  cluster_size_tbl <- tibble(
    cluster = names(cluster_sizes),
    n_cells = as.integer(cluster_sizes)
  ) %>% arrange(as.numeric(cluster))

  write_csv(cluster_size_tbl, file.path(outdir, "markers", "10_unintegrated_cluster_sizes.csv"))

  # Run markers
  markers <- FindAllMarkers(
    object = seu,
    assay = "RNA",
    slot = "data",
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    return.thresh = 0.05,
    verbose = FALSE
  )

  if (nrow(markers) == 0) {
    warning("FindAllMarkers returned zero rows. Writing empty marker files.")
    write_csv(tibble(), file.path(outdir, "markers", "10_unintegrated_markers_all_clusters.csv"))
    write_csv(tibble(), file.path(outdir, "markers", "10_unintegrated_markers_top10_per_cluster.csv"))
  } else {
    write_csv(markers, file.path(outdir, "markers", "10_unintegrated_markers_all_clusters.csv"))

    if (!"cluster" %in% colnames(markers)) {
      stop("Markers table does not contain a 'cluster' column. Columns found: ",
           paste(colnames(markers), collapse = ", "))
    }

    top10 <- markers %>%
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) %>%
      ungroup()

    write_csv(top10, file.path(outdir, "markers", "10_unintegrated_markers_top10_per_cluster.csv"))
  }

  genes_to_plot <- c(
    # endothelial
    "Pecam1","Cdh5","Kdr","Emcn","Ptprb",
    "Esm1","Apln","Angpt2","Unc5b",
    "Gja5","Efnb2","Bmx","Sox17",
    "Nr2f2","Vcam1","Vwf",
    "Cldn5","Mfsd2a","Slc2a1",
    "Mki67","Top2a","Pclaf","Cenpf",
    # photoreceptors
    "Rho","Pde6b","Sag","Gngt1","Nr2e3","Crx","Rp1",
    # neuronal
    "Syt1","Snap25","Rbfox3","Tubb3","Cadm2",
    # mural/pericyte
    "Rgs5","Pdgfrb","Cspg4","Des","Mgp","Acta2",
    # glia / Müller / astro-like
    "Rlbp1","Glul","Aqp4","Slc1a3","Sox9",
    # immune
    "Ptprc","Lyz2","Tyrobp","Ctss","C1qa",
    # RBC
    "Hbb-bs","Hba-a1","Alas2"
  )

  genes_present <- genes_to_plot[genes_to_plot %in% rownames(seu)]
  write_csv(
    tibble(gene = genes_present),
    file.path(outdir, "markers", "10_unintegrated_marker_panel_genes_present.csv")
  )

  p_umap <- DimPlot(
    seu,
    reduction = "umap.unintegrated",
    group.by = cluster_col,
    label = TRUE,
    repel = TRUE
  ) + ggtitle("Unintegrated UMAP clusters (res 0.4)")

  ggsave(
    filename = file.path("plots", "markers", "10_unintegrated_umap_clusters.png"),
    plot = p_umap,
    width = 8,
    height = 6,
    dpi = 300
  )

  batch_size <- 9
  gene_batches <- split(genes_present, ceiling(seq_along(genes_present) / batch_size))

  for (i in seq_along(gene_batches)) {
    gset <- gene_batches[[i]]
    plist <- lapply(gset, function(g) {
      FeaturePlot(
        seu,
        reduction = "umap.unintegrated",
        features = g
      ) + ggtitle(g)
    })
    p <- wrap_plots(plist, ncol = 3)
    ggsave(
      filename = file.path("plots", "markers", paste0("10_unintegrated_featureplot_batch_", i, ".png")),
      plot = p,
      width = 12,
      height = 10,
      dpi = 300
    )
  }

  cluster_tbl <- seu@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    transmute(
      cell_id = cell_id,
      sample_id = .data[[sample_col]],
      cluster = as.character(.data[[cluster_col]])
    ) %>%
    count(cluster, sample_id, name = "n") %>%
    arrange(as.numeric(cluster), sample_id)

  write_csv(cluster_tbl, file.path(outdir, "markers", "10_unintegrated_cluster_counts_by_sample.csv"))

  log_message("Stage 10_review_contaminants_unintegrated complete")
}

main()
