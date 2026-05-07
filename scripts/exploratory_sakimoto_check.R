#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_rds <- "results/objects/14_ec_unintegrated_clustered_umap.rds"

  plot_dir <- "plots/exploratory_sakimoto"
  result_dir <- "results/exploratory_sakimoto"

  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting exploratory Sakimoto-like screen")
  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  DefaultAssay(seu) <- "RNA"
  seu <- JoinLayers(seu, assay = "RNA")
  seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)

  cluster_col <- "SCT_snn_res.0.4"
  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop("Expected cluster column not found: ", cluster_col)
  }

  genes_to_check <- c("Bst1", "Procr")
  genes_present <- genes_to_check[genes_to_check %in% rownames(seu)]
  genes_missing <- setdiff(genes_to_check, genes_present)

  write_csv(
    tibble(
      gene = genes_to_check,
      present = genes_to_check %in% rownames(seu)
    ),
    file.path(result_dir, "sakimoto_marker_presence.csv")
  )

  if (length(genes_missing) > 0) {
    message_ts("Missing genes:", paste(genes_missing, collapse = ", "))
  }

  if (length(genes_present) == 0) {
    stop("Neither Bst1 nor Procr is present in the object")
  }

  avg_exp <- AverageExpression(
    seu,
    assays = "RNA",
    group.by = cluster_col,
    features = genes_present,
    slot = "data",
    verbose = FALSE
  )$RNA

  avg_df <- as.data.frame(avg_exp)
  avg_df$gene <- rownames(avg_df)
  avg_df <- avg_df %>% relocate(gene)

  write_csv(
    avg_df,
    file.path(result_dir, "sakimoto_marker_average_expression_by_cluster.csv")
  )

  total_pos_df <- FetchData(seu, vars = c(cluster_col, genes_present)) %>%
    rename(cluster = all_of(cluster_col))

  for (g in genes_present) {
    pos_summary <- total_pos_df %>%
      mutate(is_positive = .data[[g]] > 0) %>%
      group_by(cluster) %>%
      summarise(
        n_cells = n(),
        n_positive = sum(is_positive),
        pct_positive = 100 * mean(is_positive),
        avg_expression = mean(.data[[g]]),
        .groups = "drop"
      ) %>%
      arrange(as.numeric(as.character(cluster)))

    write_csv(
      pos_summary,
      file.path(result_dir, paste0("sakimoto_", g, "_positivity_by_cluster.csv"))
    )

    p_feat <- FeaturePlot(
      seu,
      features = g,
      reduction = "umap",
      raster = FALSE
    ) + ggtitle(paste0("EC-only unintegrated: ", g))

    ggsave(
      filename = file.path(plot_dir, paste0("featureplot_", g, ".png")),
      plot = p_feat,
      width = 8,
      height = 6,
      dpi = 300
    )

    p_vln <- VlnPlot(
      seu,
      features = g,
      group.by = cluster_col,
      pt.size = 0,
      raster = FALSE
    ) + NoLegend() + ggtitle(paste0(g, " by cluster"))

    ggsave(
      filename = file.path(plot_dir, paste0("vlnplot_", g, "_by_cluster.png")),
      plot = p_vln,
      width = 10,
      height = 6,
      dpi = 300
    )
  }

  if (all(c("Bst1", "Procr") %in% genes_present)) {
    coexpr_df <- FetchData(seu, vars = c(cluster_col, "orig.ident", "Bst1", "Procr")) %>%
      rename(cluster = all_of(cluster_col)) %>%
      mutate(
        Bst1_positive = Bst1 > 0,
        Procr_positive = Procr > 0,
        double_positive = Bst1_positive & Procr_positive
      )

    write_csv(
      coexpr_df,
      file.path(result_dir, "sakimoto_Bst1_Procr_cell_level_table.csv")
    )

    coexpr_summary <- coexpr_df %>%
      group_by(cluster) %>%
      summarise(
        n_cells = n(),
        n_double_positive = sum(double_positive),
        pct_double_positive = 100 * mean(double_positive),
        avg_Bst1 = mean(Bst1),
        avg_Procr = mean(Procr),
        .groups = "drop"
      ) %>%
      arrange(as.numeric(as.character(cluster)))

    write_csv(
      coexpr_summary,
      file.path(result_dir, "sakimoto_Bst1_Procr_double_positive_by_cluster.csv")
    )

    p_blend <- FeaturePlot(
      seu,
      features = c("Bst1", "Procr"),
      reduction = "umap",
      blend = TRUE,
      raster = FALSE
    ) + ggtitle("EC-only unintegrated: Bst1 / Procr blend")

    ggsave(
      filename = file.path(plot_dir, "featureplot_Bst1_Procr_blend.png"),
      plot = p_blend,
      width = 8,
      height = 6,
      dpi = 300
    )
  }

  message_ts("Exploratory Sakimoto-like screen complete")
}

main()
