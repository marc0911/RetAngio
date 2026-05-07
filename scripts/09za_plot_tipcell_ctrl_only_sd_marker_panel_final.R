suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

find_first_existing <- function(paths) {
  paths <- unique(paths)
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

find_cluster_col <- function(meta) {
  candidates <- c(
    "tip_ctrl_res_0_20",
    "seurat_clusters",
    "tip_clusters",
    "local_cluster",
    "cluster"
  )
  hit <- candidates[candidates %in% colnames(meta)]
  if (length(hit) == 0) {
    stop(
      "Could not find a local cluster column. Available columns: ",
      paste(colnames(meta), collapse = ", ")
    )
  }
  hit[1]
}

main <- function() {
  message_ts("Starting stage 09za_plot_tipcell_ctrl_only_sd_marker_panel_final")

  candidate_objects <- c(
    "results/objects/09o_rpca_tipcell_ctrl_only.rds",
    "results/objects/09p_rpca_tipcell_ctrl_only.rds",
    list.files("results/objects", pattern = "tip.*ctrl.*\\.rds$", full.names = TRUE)
  )

  input_rds <- find_first_existing(candidate_objects)
  if (is.na(input_rds)) {
    stop("No suitable control-only tip-cell object found in results/objects/")
  }

  message_ts("Loading object:", input_rds)
  seu <- readRDS(input_rds)

  output_dir <- "results/review"
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  assay_to_use <- if ("RNA" %in% Assays(seu)) "RNA" else DefaultAssay(seu)
  DefaultAssay(seu) <- assay_to_use
  message_ts("Using assay:", assay_to_use)

  if (assay_to_use == "RNA") {
    message_ts("Attempting to join RNA layers if needed")
    try({
      seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
    }, silent = TRUE)
  }

  cluster_col <- find_cluster_col(seu@meta.data)
  message_ts("Using cluster column:", cluster_col)

  cluster_values <- seu@meta.data[[cluster_col]]
  if (is.numeric(cluster_values) || all(grepl("^[0-9]+$", as.character(cluster_values)))) {
    cluster_levels <- as.character(sort(unique(as.numeric(as.character(cluster_values)))))
  } else {
    cluster_levels <- unique(as.character(cluster_values))
  }

  seu[[cluster_col]] <- factor(as.character(cluster_values), levels = cluster_levels)
  Idents(seu) <- seu[[cluster_col]][, 1]

  marker_groups <- list(
    "S-tip" = c("Esm1", "Igfbp3", "Apln", "Angpt2", "Col15a1"),
    "D-tip" = c("Spock2", "Apod", "Cldn5", "Slc38a5", "Itm2a")
  )

  genes_all <- unlist(marker_groups, use.names = FALSE)
  genes_present <- genes_all[genes_all %in% rownames(seu)]

  genes_checked <- data.frame(
    gene = genes_all,
    present = genes_all %in% rownames(seu),
    stringsAsFactors = FALSE
  )

  write.csv(
    genes_checked,
    file = file.path(output_dir, "09za_rpca_tip_ctrl_sd_marker_panel_genes_checked.csv"),
    row.names = FALSE
  )

  marker_groups_present <- lapply(marker_groups, function(x) x[x %in% rownames(seu)])

  message_ts(
    "Genes present:",
    paste(unlist(marker_groups_present, use.names = FALSE), collapse = ", ")
  )

  p <- DotPlot(
    object = seu,
    features = marker_groups_present,
    cols = c("#2166AC", "white", "#B2182B"),
    dot.scale = 10,
    scale = TRUE
  ) +
    labs(
      title = NULL,
      x = "Genes",
      y = "Local cluster"
    ) +
    theme_classic(base_size = 18) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title.x = element_text(face = "bold", size = 18),
      axis.title.y = element_text(face = "bold", size = 18),
      axis.text = element_text(face = "bold", size = 14),
      strip.text = element_text(face = "bold", size = 16),
      strip.background = element_rect(fill = "grey90", colour = "black"),
      plot.title = element_blank(),
      legend.title = element_text(face = "bold", size = 16),
      legend.text = element_text(size = 13)
    )

  ggsave(
    filename = file.path(output_dir, "09za_rpca_tip_ctrl_sd_marker_panel_final.pdf"),
    plot = p,
    width = 14,
    height = 9
  )

  ggsave(
    filename = file.path(output_dir, "09za_rpca_tip_ctrl_sd_marker_panel_final.png"),
    plot = p,
    width = 14,
    height = 9,
    dpi = 300
  )

  message_ts("Stage 09za_plot_tipcell_ctrl_only_sd_marker_panel_final complete")
}

main()
