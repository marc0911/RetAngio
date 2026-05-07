suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Cairo)
})

log_message <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_rds <- "results/objects/09h_ec_after_soupx_rpca_annotated.rds"
  output_dir <- "plots/09r_ec_after_soupx_rpca_curated_panel"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  log_message("Starting stage 09r_plot_curated_overall_ec_panel_rpca")
  log_message("Loading object:", input_rds)
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

  annotation_col <- "annotation_my_label"
  log_message("Using annotation column:", annotation_col)

  label_order <- c(
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

  seu@meta.data[[annotation_col]] <- factor(
    seu@meta.data[[annotation_col]],
    levels = label_order
  )

  assay_to_use <- if ("RNA" %in% names(seu@assays)) "RNA" else "SCT"
  DefaultAssay(seu) <- assay_to_use
  log_message("Using assay:", assay_to_use)

  if (assay_to_use == "RNA") {
    assay_obj <- seu[["RNA"]]

    if (inherits(assay_obj, "Assay5")) {
      lyr <- Layers(assay_obj)

      if (length(lyr) > 1) {
        log_message("Joining RNA layers")
        seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
        assay_obj <- seu[["RNA"]]
        lyr <- Layers(assay_obj)
      }

      if (!("data" %in% lyr)) {
        log_message("Normalizing RNA assay")
        seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)
      }
    } else {
      if (nrow(assay_obj@data) == 0) {
        log_message("Normalizing RNA assay")
        seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)
      }
    }
  }

  gene_panel <- c(
    "Bmx", "Gja5", "Efnb2",
    "Cldn5", "Mfsd2a", "Slc2a1", "Flt1",
    "Slc38a5", "Vcam1", "Nr2f2",
    "Apln", "Angpt2", "Dll4", "Cxcr4", "Esm1",
    "Mki67", "Top2a", "Cdk1",
    "Ifitm3", "Isg15", "Bst2",
    "Aqp1", "Six3", "Slco1a4",
    "Cfh", "Igfbp4", "Fn1"
  )

  genes_present <- gene_panel[gene_panel %in% rownames(seu)]
  genes_missing <- setdiff(gene_panel, genes_present)

  write.csv(
    data.frame(
      gene = gene_panel,
      present_in_object = gene_panel %in% genes_present
    ),
    file = file.path(output_dir, "09r_rpca_curated_panel_genes_checked.csv"),
    row.names = FALSE
  )

  if (length(genes_missing) > 0) {
    log_message("Missing genes:", paste(genes_missing, collapse = ", "))
  }
  log_message("Genes present:", paste(genes_present, collapse = ", "))

  p <- DotPlot(
    object = seu,
    features = genes_present,
    group.by = annotation_col,
    dot.scale = 6
  ) +
    scale_y_discrete(limits = rev(label_order)) +
    scale_colour_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "Scaled average\nexpression"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(face = "bold"),
      panel.grid.major = element_line(linewidth = 0.2),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(
      title = "RPCA overall EC object - curated marker panel",
      x = NULL,
      y = NULL
    )

  p_nolegend <- p + theme(legend.position = "none")

  CairoPDF(file.path(output_dir, "09r_rpca_overall_ec_curated_dotplot.pdf"), width = 16, height = 8)
  print(p)
  dev.off()

  CairoPNG(file.path(output_dir, "09r_rpca_overall_ec_curated_dotplot.png"), width = 16, height = 8, units = "in", dpi = 300)
  print(p)
  dev.off()

  CairoPDF(file.path(output_dir, "09r_rpca_overall_ec_curated_dotplot_NOlegend.pdf"), width = 16, height = 8)
  print(p_nolegend)
  dev.off()

  CairoPNG(file.path(output_dir, "09r_rpca_overall_ec_curated_dotplot_NOlegend.png"), width = 16, height = 8, units = "in", dpi = 300)
  print(p_nolegend)
  dev.off()

  log_message("Saved plots to:", output_dir)
  log_message("Stage 09r_plot_curated_overall_ec_panel_rpca complete")
}

main()
