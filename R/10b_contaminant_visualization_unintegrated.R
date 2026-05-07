source("R/00_utils.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

main <- function() {
  log_message("Starting stage 10b_contaminant_visualization_unintegrated")

  in_rds <- "results/objects/08_unintegrated_clustered_umap.rds"
  outdir <- "plots/contaminants"

  ensure_dir("plots")
  ensure_dir(outdir)

  seu <- readRDS(in_rds)
  log_message("Loaded object: ", in_rds)

  cluster_col <- "unintegrated_snn_res.0.4"
  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop("Missing cluster column: ", cluster_col)
  }

  DefaultAssay(seu) <- "RNA"
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)

  Idents(seu) <- seu@meta.data[[cluster_col]]

  # Base UMAP
  p_umap <- DimPlot(
    seu,
    reduction = "umap.unintegrated",
    group.by = cluster_col,
    label = TRUE,
    repel = TRUE
  ) + ggtitle("Unintegrated UMAP clusters (res 0.4)")

  ggsave(
    file.path(outdir, "10b_umap_clusters.png"),
    p_umap, width = 8, height = 6, dpi = 300
  )

  # Contaminant marker groups
  marker_groups <- list(
    photoreceptor = c("Rho", "Pde6b", "Sag", "Gngt1", "Gnat1", "Nr2e3", "Crx"),
    neuronal = c("Syt1", "Snap25", "Rbfox3", "Tubb3", "Cadm2", "Sncg"),
    mural_pericyte = c("Cspg4", "Pdgfrb", "Rgs5", "Ednra", "Kcnj8", "Notch3"),
    glia_muller = c("Rlbp1", "Glul", "Aqp4", "Slc1a3", "Sox9"),
    immune = c("Ptprc", "Lyz2", "Tyrobp", "Ctss", "C1qa"),
    rbc = c("Hbb-bs", "Hba-a1", "Alas2"),
    endothelial = c("Pecam1", "Cdh5", "Kdr", "Ptprb", "Cldn5", "Kdr")
  )

  # FeaturePlots by group
  for (nm in names(marker_groups)) {
    genes <- marker_groups[[nm]]
    genes <- genes[genes %in% rownames(seu)]
    if (length(genes) == 0) next

    plist <- lapply(genes, function(g) {
      FeaturePlot(
        seu,
        reduction = "umap.unintegrated",
        features = g
      ) + ggtitle(g)
    })

    p <- wrap_plots(plist, ncol = 3)
    ggsave(
      file.path(outdir, paste0("10b_featureplot_", nm, ".png")),
      p, width = 12, height = 8, dpi = 300
    )
  }

  # DotPlot summary for decision making
  dot_genes <- unique(c(
    marker_groups$endothelial,
    marker_groups$photoreceptor,
    marker_groups$neuronal,
    marker_groups$mural_pericyte,
    marker_groups$glia_muller,
    marker_groups$immune,
    marker_groups$rbc
  ))
  dot_genes <- dot_genes[dot_genes %in% rownames(seu)]

  p_dot <- DotPlot(
    seu,
    features = dot_genes,
    assay = "RNA"
  ) +
    RotatedAxis() +
    ggtitle("Contaminant and endothelial markers by cluster")

  ggsave(
    file.path(outdir, "10b_dotplot_contaminants_vs_endothelial.png"),
    p_dot, width = 16, height = 8, dpi = 300
  )

  log_message("Stage 10b_contaminant_visualization_unintegrated complete")
}

main()
