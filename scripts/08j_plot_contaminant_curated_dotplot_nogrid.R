#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(Cairo)
})

log_message <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", paste0(..., collapse = ""), "\n")
}

main <- function() {
  stage_name <- "08j_plot_contaminant_curated_dotplot_nogrid"

  input_rds <- "results/objects/08_unintegrated_clustered_umap.rds"
  output_dir <- "plots/08j_precontaminant_curated_contaminant_dotplot_nogrid"
  review_dir <- "results/review"
  cluster_col <- "seurat_clusters"

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(review_dir, recursive = TRUE, showWarnings = FALSE)

  log_message("Starting stage ", stage_name)
  log_message("Loading object: ", input_rds)
  seu <- readRDS(input_rds)

  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop(
      "Cluster column not found: ", cluster_col,
      "\nAvailable metadata columns: ",
      paste(colnames(seu@meta.data), collapse = ", ")
    )
  }

  assay_to_use <- if ("RNA" %in% Assays(seu)) "RNA" else DefaultAssay(seu)
  log_message("Using assay: ", assay_to_use)
  DefaultAssay(seu) <- assay_to_use

  if (assay_to_use == "RNA") {
    assay_obj <- seu[["RNA"]]
    if (inherits(assay_obj, "Assay5")) {
      lyr <- Layers(assay_obj)
      if (length(lyr) > 1) {
        log_message("Joining RNA layers")
        seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
      }
      if (!("data" %in% Layers(seu[["RNA"]]))) {
        log_message("Normalizing RNA assay")
        seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)
      }
    } else {
      if (nrow(seu[["RNA"]]@data) == 0) {
        log_message("Normalizing RNA assay")
        seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)
      }
    }
  }

  cluster_values <- as.character(seu@meta.data[[cluster_col]])
  cluster_levels <- unique(cluster_values)

  if (all(grepl("^[0-9]+$", cluster_levels))) {
    cluster_levels <- as.character(sort(as.numeric(cluster_levels)))
  } else {
    cluster_levels <- sort(cluster_levels)
  }

  seu$cluster_plot_08j <- factor(cluster_values, levels = cluster_levels)
  Idents(seu) <- "cluster_plot_08j"

  log_message("Using cluster column: ", cluster_col)
  log_message("Cluster levels: ", paste(cluster_levels, collapse = ", "))

  marker_panel <- list(
    "endothelial" = c("Pecam1", "Cdh5", "Kdr"),
    "immune / myeloid" = c("Ptprc", "Lyz2", "C1qa"),
    "photoreceptor" = c("Rho", "Rcvrn", "Pde6b"),
    "neuronal / retina" = c("Snap25", "Tubb3", "Stmn2"),
    "glia / Muller" = c("Rlbp1", "Slc1a3", "Gfap"),
    "mural / pericyte" = c("Pdgfrb", "Rgs5", "Kcnj8"),
    "RBC" = c("Hbb-bs", "Hba-a1", "Alas2")
  )

  all_genes <- unique(unlist(marker_panel))
  genes_present <- all_genes[all_genes %in% rownames(seu)]
  genes_missing <- setdiff(all_genes, genes_present)

  log_message("Genes present: ", paste(genes_present, collapse = ", "))
  if (length(genes_missing) > 0) {
    log_message("Genes missing: ", paste(genes_missing, collapse = ", "))
  }

  write_csv(
    data.frame(
      gene = all_genes,
      present_in_object = all_genes %in% rownames(seu)
    ),
    file.path(review_dir, "08j_precontaminant_curated_contaminant_dotplot_genes_checked.csv")
  )

  marker_panel_present <- lapply(marker_panel, function(x) x[x %in% rownames(seu)])
  marker_panel_present <- marker_panel_present[sapply(marker_panel_present, length) > 0]

  dp <- DotPlot(
    object = seu,
    features = marker_panel_present,
    assay = assay_to_use,
    dot.scale = 8
  )

  dp_data <- dp$data
  gene_order <- unlist(marker_panel_present, use.names = FALSE)
  dp_data$features.plot <- factor(dp_data$features.plot, levels = gene_order)
  dp_data$id <- factor(dp_data$id, levels = rev(cluster_levels))

  write_csv(
    dp_data,
    file.path(review_dir, "08j_precontaminant_curated_contaminant_dotplot_data.csv")
  )

  p <- ggplot(dp_data, aes(x = features.plot, y = id)) +
    geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
    facet_grid(
      . ~ feature.groups,
      scales = "free_x",
      space = "free_x",
      switch = "x"
    ) +
    scale_size(
      name = "Percent Expressed",
      range = c(0, 8),
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100)
    ) +
    scale_color_gradient2(
      name = "Average Expression",
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0
    ) +
    labs(
      title = "Pre-contaminant-removal object - curated contaminant marker panel",
      x = "Genes",
      y = "Cluster"
    ) +
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
      axis.title = element_text(face = "bold", size = 18),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 13, color = "black"),
      axis.text.y = element_text(face = "bold", size = 13, color = "black"),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", linewidth = 0.8),
      strip.background = element_rect(fill = "white", color = "black", linewidth = 1),
      strip.text = element_text(face = "bold", size = 14),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12)
    )

  p_nolegend <- p + theme(legend.position = "none")

  CairoPDF(
    file.path(output_dir, "08j_precontaminant_curated_contaminant_dotplot_nogrid.pdf"),
    width = 22,
    height = 10
  )
  print(p)
  dev.off()

  CairoPNG(
    file.path(output_dir, "08j_precontaminant_curated_contaminant_dotplot_nogrid.png"),
    width = 22,
    height = 10,
    units = "in",
    dpi = 300
  )
  print(p)
  dev.off()

  CairoPDF(
    file.path(output_dir, "08j_precontaminant_curated_contaminant_dotplot_nogrid_NOlegend.pdf"),
    width = 22,
    height = 10
  )
  print(p_nolegend)
  dev.off()

  CairoPNG(
    file.path(output_dir, "08j_precontaminant_curated_contaminant_dotplot_nogrid_NOlegend.png"),
    width = 22,
    height = 10,
    units = "in",
    dpi = 300
  )
  print(p_nolegend)
  dev.off()

  log_message("Stage ", stage_name, " complete")
}

main()
