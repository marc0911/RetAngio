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

safe_join_and_normalize_rna <- function(seu) {
  if (!"RNA" %in% Assays(seu)) return(seu)

  DefaultAssay(seu) <- "RNA"
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

  seu
}

find_input_rds <- function() {
  candidates <- c(
    "results/objects/09o_rpca_tipcell_ctrl_only_resolution_sweep.rds",
    "results/objects/09o_rpca_tipcell_ctrl_only_res_0_20.rds",
    "results/objects/09p_rpca_tipcell_ctrl_only_resolution_sweep.rds",
    "results/objects/09p_rpca_tipcell_ctrl_only_res_0_20.rds",
    "results/objects/09q_rpca_tipcell_ctrl_only_res_0_20.rds"
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop(
      "Could not find ctrl-only tip-cell object.\nChecked:\n",
      paste(candidates, collapse = "\n")
    )
  }
  existing[1]
}

find_cluster_col <- function(meta) {
  preferred <- c(
    "tip_ctrl_res_0_20",
    "tipcell_ctrl_res_0_20",
    "ctrl_tip_res_0_20",
    "tip_res_0_20",
    "seurat_clusters"
  )
  hit <- preferred[preferred %in% colnames(meta)]
  if (length(hit) == 0) {
    stop(
      "Could not find a suitable cluster column.\nAvailable columns:\n",
      paste(colnames(meta), collapse = ", ")
    )
  }
  hit[1]
}

get_cluster_levels <- function(x) {
  lev <- unique(as.character(x))
  if (all(grepl("^[0-9]+$", lev))) {
    lev <- as.character(sort(as.numeric(lev)))
  } else {
    lev <- sort(lev)
  }
  lev
}

main <- function() {
  stage_name <- "09zc_fix_ctrl_tip_sd_dotplot_labels_only"
  output_dir <- "plots/09zc_fix_ctrl_tip_sd_dotplot_labels_only"
  review_dir <- "results/review"

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(review_dir, recursive = TRUE, showWarnings = FALSE)

  input_rds <- find_input_rds()
  log_message("Starting stage ", stage_name)
  log_message("Loading object: ", input_rds)

  seu <- readRDS(input_rds)

  cluster_col <- find_cluster_col(seu@meta.data)
  log_message("Using cluster column: ", cluster_col)

  assay_to_use <- if ("RNA" %in% Assays(seu)) "RNA" else DefaultAssay(seu)
  log_message("Using assay: ", assay_to_use)
  DefaultAssay(seu) <- assay_to_use

  if (assay_to_use == "RNA") {
    seu <- safe_join_and_normalize_rna(seu)
  }

  cluster_values <- as.character(seu@meta.data[[cluster_col]])
  cluster_levels <- get_cluster_levels(cluster_values)
  log_message("Cluster levels: ", paste(cluster_levels, collapse = ", "))

  seu$ctrl_tip_plot_clusters_09zc <- factor(cluster_values, levels = cluster_levels)
  Idents(seu) <- "ctrl_tip_plot_clusters_09zc"

  marker_groups <- list(
    "S-tip" = c("Esm1", "Igfbp3", "Apln", "Angpt2", "Col15a1"),
    "D-tip" = c("Spock2", "Apod", "Cldn5", "Slc38a5", "Itm2a")
  )

  all_genes <- unique(unlist(marker_groups, use.names = FALSE))
  genes_present <- all_genes[all_genes %in% rownames(seu)]
  genes_missing <- setdiff(all_genes, genes_present)

  log_message("Genes present: ", paste(genes_present, collapse = ", "))
  if (length(genes_missing) > 0) {
    log_message("Genes missing: ", paste(genes_missing, collapse = ", "))
  }

  write_csv(
    data.frame(
      group = rep(names(marker_groups), lengths(marker_groups)),
      gene = unlist(marker_groups, use.names = FALSE),
      present_in_object = unlist(marker_groups, use.names = FALSE) %in% rownames(seu)
    ),
    file.path(review_dir, "09zc_ctrl_tip_SD_marker_panel_used.csv")
  )

  marker_groups_present <- lapply(marker_groups, function(x) x[x %in% rownames(seu)])
  marker_groups_present <- marker_groups_present[sapply(marker_groups_present, length) > 0]

  dp <- DotPlot(
    object = seu,
    features = marker_groups_present,
    assay = assay_to_use,
    dot.scale = 8
  )

  dp_data <- dp$data
  gene_order <- unlist(marker_groups_present, use.names = FALSE)
  dp_data$features.plot <- factor(dp_data$features.plot, levels = gene_order)
  dp_data$id <- factor(dp_data$id, levels = rev(cluster_levels))

  write_csv(
    dp_data,
    file.path(review_dir, "09zc_ctrl_tip_SD_dotplot_data.csv")
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
      title = NULL,
      x = "Genes",
      y = "Local cluster"
    ) +
    theme_classic(base_size = 18) +
    theme(
      plot.title = element_blank(),
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

  CairoPDF(
    file.path(output_dir, "09zc_ctrl_tip_SD_marker_panel.pdf"),
    width = 13,
    height = 8
  )
  print(p)
  dev.off()

  CairoPNG(
    file.path(output_dir, "09zc_ctrl_tip_SD_marker_panel.png"),
    width = 13,
    height = 8,
    units = "in",
    dpi = 300
  )
  print(p)
  dev.off()

  log_message("Stage ", stage_name, " complete")
}

main()
