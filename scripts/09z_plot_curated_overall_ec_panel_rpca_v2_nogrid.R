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
  stage_name <- "09z_plot_curated_overall_ec_panel_rpca_v2_nogrid"

  input_rds <- "results/objects/09h_ec_after_soupx_rpca_annotated.rds"
  output_dir <- "plots/09z_rpca_overall_ec_curated_panel_v2_nogrid"
  review_dir <- "results/review"

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(review_dir, recursive = TRUE, showWarnings = FALSE)

  log_message("Starting stage ", stage_name)
  log_message("Loading object: ", input_rds)
  seu <- readRDS(input_rds)

  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    stop("Missing seurat_clusters in metadata.")
  }

  rpca_label_map <- c(
    "0"  = "Venous x capillary",
    "1"  = "Capillary 1",
    "2"  = "Capillary 2",
    "3"  = "Arterial x capillary",
    "4"  = "Proliferative 1",
    "5"  = "Tip cell",
    "6"  = "Arterial",
    "7"  = "EC1",
    "8"  = "Proliferative 2",
    "9"  = "Uncertain",
    "10" = "Proliferative 3",
    "11" = "Inflammatory EC",
    "12" = "EC2"
  )

  seu$annotation_my_label <- unname(rpca_label_map[as.character(seu$seurat_clusters)])

  if (any(is.na(seu$annotation_my_label))) {
    missing_clusters <- sort(unique(as.character(seu$seurat_clusters[is.na(seu$annotation_my_label)])))
    stop("Some clusters are not mapped: ", paste(missing_clusters, collapse = ", "))
  }

  annotation_order <- c(
    "Arterial",
    "Arterial x capillary",
    "Capillary 1",
    "Capillary 2",
    "Venous x capillary",
    "Tip cell",
    "Proliferative 1",
    "Proliferative 2",
    "Proliferative 3",
    "Inflammatory EC",
    "EC1",
    "EC2",
    "Uncertain"
  )

  seu$annotation_my_label <- factor(seu$annotation_my_label, levels = annotation_order)
  Idents(seu) <- "annotation_my_label"

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

  marker_panel <- list(
    "arterial" = c("Bmx", "Gja5", "Efnb2"),
    "capillary/BRB" = c("Cldn5", "Mfsd2a", "Slc2a1"),
    "venous x capillary" = c("Apod", "Colec12", "Adm"),
    "tip" = c("Apln", "Angpt2", "Esm1"),
    "proliferative" = c("Mki67", "Top2a"),
    "inflammatory" = c("Isg15", "Bst2"),
    "EC1" = c("Aqp1", "Six3"),
    "EC2" = c("Ntf3", "Ihh", "Rarb")
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
    file.path(review_dir, "09z_rpca_curated_panel_v2_nogrid_genes_checked.csv")
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
  dp_data$id <- factor(dp_data$id, levels = rev(annotation_order))

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
      title = "RPCA overall EC object - curated marker panel (v2)",
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

  CairoPDF(file.path(output_dir, "09z_rpca_overall_ec_curated_dotplot_v2_nogrid.pdf"), width = 20, height = 10)
  print(p)
  dev.off()

  CairoPNG(
    file.path(output_dir, "09z_rpca_overall_ec_curated_dotplot_v2_nogrid.png"),
    width = 20,
    height = 10,
    units = "in",
    dpi = 300
  )
  print(p)
  dev.off()

  CairoPDF(file.path(output_dir, "09z_rpca_overall_ec_curated_dotplot_v2_nogrid_NOlegend.pdf"), width = 20, height = 10)
  print(p_nolegend)
  dev.off()

  CairoPNG(
    file.path(output_dir, "09z_rpca_overall_ec_curated_dotplot_v2_nogrid_NOlegend.png"),
    width = 20,
    height = 10,
    units = "in",
    dpi = 300
  )
  print(p_nolegend)
  dev.off()

  log_message("Stage ", stage_name, " complete")
}

main()
