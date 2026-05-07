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

get_cluster_levels <- function(x) {
  lev <- unique(as.character(x))
  if (all(grepl("^[0-9]+$", lev))) {
    lev <- as.character(sort(as.numeric(lev)))
  } else {
    lev <- sort(lev)
  }
  lev
}

make_dotplot <- function(seu, marker_groups, cluster_levels, assay_to_use,
                         output_prefix, output_dir, review_dir,
                         width, height) {
  all_genes <- unlist(marker_groups, use.names = FALSE)

  marker_table <- data.frame(
    group = rep(names(marker_groups), lengths(marker_groups)),
    gene = all_genes,
    present_in_object = all_genes %in% rownames(seu),
    stringsAsFactors = FALSE
  )

  write_csv(
    marker_table,
    file.path(review_dir, paste0(output_prefix, "_marker_panel_used.csv"))
  )

  genes_missing <- marker_table$gene[!marker_table$present_in_object]
  if (length(genes_missing) > 0) {
    log_message("Genes missing for ", output_prefix, ": ", paste(genes_missing, collapse = ", "))
  }

  marker_groups_present <- lapply(marker_groups, function(x) x[x %in% rownames(seu)])
  marker_groups_present <- marker_groups_present[sapply(marker_groups_present, length) > 0]

  log_message(
    "Genes present for ", output_prefix, ": ",
    paste(unlist(marker_groups_present, use.names = FALSE), collapse = ", ")
  )

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
    file.path(review_dir, paste0(output_prefix, "_dotplot_data.csv"))
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

  p_nolegend <- p + theme(legend.position = "none")

  CairoPDF(file.path(output_dir, paste0(output_prefix, ".pdf")), width = width, height = height)
  print(p)
  dev.off()

  CairoPNG(
    file.path(output_dir, paste0(output_prefix, ".png")),
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )
  print(p)
  dev.off()

  CairoPDF(file.path(output_dir, paste0(output_prefix, "_NOlegend.pdf")), width = width, height = height)
  print(p_nolegend)
  dev.off()

  CairoPNG(
    file.path(output_dir, paste0(output_prefix, "_NOlegend.png")),
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )
  print(p_nolegend)
  dev.off()
}

main <- function() {
  stage_name <- "09zd_fix_overall_tip_sd_pathological_dotplots"

  input_rds <- "results/objects/09k_rpca_tipcell_resolution_sweep.rds"
  cluster_col <- "tip_res_0_20"
  output_dir <- "plots/09zd_fix_overall_tip_sd_pathological_dotplots"
  review_dir <- "results/review"

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
    seu <- safe_join_and_normalize_rna(seu)
  }

  cluster_values <- as.character(seu@meta.data[[cluster_col]])
  cluster_levels <- get_cluster_levels(cluster_values)

  seu$overall_tip_plot_clusters_09zd <- factor(cluster_values, levels = cluster_levels)
  Idents(seu) <- "overall_tip_plot_clusters_09zd"

  log_message("Using cluster column: ", cluster_col)
  log_message("Cluster levels: ", paste(cluster_levels, collapse = ", "))

  sd_marker_groups <- list(
    "S-tip" = c("Esm1", "Igfbp3", "Apln", "Angpt2", "Col15a1"),
    "D-tip" = c("Spock2", "Apod", "Cldn5", "Slc38a5", "Itm2a")
  )

  sd_path_marker_groups <- list(
    "S-tip" = c("Esm1", "Igfbp3", "Apln", "Angpt2", "Col15a1"),
    "D-tip" = c("Spock2", "Apod", "Cldn5", "Slc38a5", "Itm2a"),
    "Pathological 1" = c("Pgf", "Bnip3", "Fabp5"),
    "Pathological 2" = c("Vegfa", "Bmp6", "Ncam1"),
    "Pathological 3" = c("Plvap", "Ednrb", "Ptn")
  )

  make_dotplot(
    seu = seu,
    marker_groups = sd_marker_groups,
    cluster_levels = cluster_levels,
    assay_to_use = assay_to_use,
    output_prefix = "09zd_overall_tip_SD_marker_panel",
    output_dir = output_dir,
    review_dir = review_dir,
    width = 13,
    height = 8
  )

  make_dotplot(
    seu = seu,
    marker_groups = sd_path_marker_groups,
    cluster_levels = cluster_levels,
    assay_to_use = assay_to_use,
    output_prefix = "09zd_overall_tip_SD_plus_curated_pathological_marker_panel",
    output_dir = output_dir,
    review_dir = review_dir,
    width = 20,
    height = 8
  )

  log_message("Stage ", stage_name, " complete")
}

main()
