suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(Cairo)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_rds <- Sys.getenv("INPUT_RDS", unset = "results/objects/09k_rpca_tipcell_resolution_sweep.rds")
  output_dir <- "results/review"
  cluster_col <- "tip_res_0_20"

  gene_panels <- list(
    "S-tip-like" = c("Esm1", "Apln", "Angpt2"),
    "D-side / D-tip-like" = c("Cldn5", "Igf1r", "Pdgfb")
  )

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 09w_plot_tipcell_overall_curated_panel_unscaled")
  message_ts("Loading object:", input_rds)

  seu <- readRDS(input_rds)

  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop(
      "Required cluster column not found: ", cluster_col,
      "\nAvailable metadata columns: ",
      paste(colnames(seu@meta.data), collapse = ", ")
    )
  }

  assay_to_use <- if ("RNA" %in% Assays(seu)) "RNA" else DefaultAssay(seu)
  DefaultAssay(seu) <- assay_to_use
  message_ts("Using assay:", assay_to_use)

  if (assay_to_use == "RNA") {
    assay_obj <- seu[["RNA"]]
    if (inherits(assay_obj, "Assay5")) {
      lyr <- Layers(assay_obj)
      if (length(lyr) > 1) {
        message_ts("Joining RNA layers")
        seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
      }
      if (!("data" %in% Layers(seu[["RNA"]]))) {
        message_ts("Normalizing RNA assay")
        seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)
      }
    } else {
      if (nrow(seu[["RNA"]]@data) == 0) {
        message_ts("Normalizing RNA assay")
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

  seu$tip_cluster_plot_09w <- factor(cluster_values, levels = cluster_levels)
  Idents(seu) <- "tip_cluster_plot_09w"

  message_ts("Using cluster column:", cluster_col)
  message_ts("Cluster levels:", paste(cluster_levels, collapse = ", "))

  genes_requested <- unlist(gene_panels, use.names = FALSE)
  genes_present <- genes_requested[genes_requested %in% rownames(seu)]
  genes_missing <- setdiff(genes_requested, genes_present)

  write.csv(
    data.frame(
      gene = genes_requested,
      present = genes_requested %in% rownames(seu),
      stringsAsFactors = FALSE
    ),
    file = file.path(output_dir, "09w_rpca_tip_overall_unscaled_genes_checked.csv"),
    row.names = FALSE
  )

  if (length(genes_missing) > 0) {
    stop("Missing genes: ", paste(genes_missing, collapse = ", "))
  }

  message_ts("Genes present:", paste(genes_present, collapse = ", "))

  p <- DotPlot(
    object = seu,
    features = gene_panels,
    assay = assay_to_use,
    scale = FALSE,
    dot.scale = 10
  )

  p$data$id <- factor(as.character(p$data$id), levels = rev(cluster_levels))

  write.csv(
    p$data,
    file = file.path(output_dir, "09w_rpca_tip_overall_unscaled_dotplot_data.csv"),
    row.names = FALSE
  )

  p <- p +
    scale_color_gradient(
      low = "grey95",
      high = "firebrick3",
      name = "Average Expression"
    ) +
    labs(
      title = "RPCA overall tip-cell object - curated marker panel (unscaled)",
      x = "Genes",
      y = "Local cluster",
      size = "Percent Expressed"
    ) +
    theme_classic(base_size = 18) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = element_text(face = "bold", color = "black"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 22),
      strip.background = element_rect(fill = "grey95", colour = "black"),
      strip.text = element_text(face = "bold", size = 16),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 12),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8)
    )

  ggsave(
    filename = file.path(output_dir, "09w_rpca_tip_overall_curated_dotplot_unscaled.pdf"),
    plot = p,
    width = 12,
    height = 7,
    units = "in",
    device = cairo_pdf
  )

  CairoPNG(
    file.path(output_dir, "09w_rpca_tip_overall_curated_dotplot_unscaled.png"),
    width = 12,
    height = 7,
    units = "in",
    dpi = 300
  )
  print(p)
  dev.off()

  message_ts("Stage 09w_plot_tipcell_overall_curated_panel_unscaled complete")
}

main()
