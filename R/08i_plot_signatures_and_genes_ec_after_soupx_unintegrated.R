#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(readxl)
  library(ggplot2)
  library(tidyr)
  library(Matrix)
})

source("R/plotting_helpers.R")

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

make_safe_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

clean_gene_vector <- function(x) {
  x <- unlist(strsplit(x, ","))
  x <- trimws(x)
  x <- x[nchar(x) > 0]
  unique(x)
}

main <- function() {
  input_rds <- "results/objects/08h_ec_after_soupx_unintegrated_signature_scored.rds"
  input_sig_used <- "results/signatures/08h_ec_after_soupx_unintegrated_signatures_used.csv"
  input_cluster_scores <- "results/signatures/08h_ec_after_soupx_unintegrated_cluster_mean_scores.csv"
  input_xlsx <- "reference/RetAngio_annotation_master_v2.xlsx"

  outdir <- "plots/08i_ec_after_soupx_unintegrated_signatures_and_genes"

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 08i_plot_signatures_and_genes_ec_after_soupx_unintegrated")

  seu <- readRDS(input_rds)
  sig_used <- read_csv(input_sig_used, show_col_types = FALSE)
  cluster_scores <- read_csv(input_cluster_scores, show_col_types = FALSE)
  sig_tbl <- readxl::read_excel(input_xlsx, sheet = "category_signatures")

  message_ts("Loaded object:", input_rds)
  message_ts("Loaded signatures used:", input_sig_used)
  message_ts("Loaded cluster mean scores:", input_cluster_scores)
  message_ts("Loaded workbook sheet: category_signatures")

  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    seu$seurat_clusters <- as.character(Idents(seu))
  }
  Idents(seu) <- "seurat_clusters"

  reduction_to_use <- if ("umap.unintegrated" %in% Reductions(seu)) {
    "umap.unintegrated"
  } else if ("umap" %in% Reductions(seu)) {
    "umap"
  } else {
    stop("No UMAP reduction found in object.")
  }

  sig_plot_tbl <- sig_used %>%
    filter(signature_set == "broader") %>%
    mutate(category_clean = make_safe_name(category))

  score_cols <- sig_plot_tbl$score_name
  missing_scores <- setdiff(score_cols, colnames(seu@meta.data))
  if (length(missing_scores) > 0) {
    stop(paste("Missing score columns in object:", paste(missing_scores, collapse = ", ")))
  }

  if (!"RNA" %in% Assays(seu)) {
    stop("RNA assay not found in object.")
  }

  DefaultAssay(seu) <- "RNA"
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  seu <- NormalizeData(seu, verbose = FALSE)

  gene_tbl <- sig_tbl %>%
    filter(category %in% sig_plot_tbl$category) %>%
    transmute(category, minimal = recommended_minimal_signature)

  main_gene_map <- lapply(seq_len(nrow(gene_tbl)), function(i) {
    genes <- clean_gene_vector(as.character(gene_tbl$minimal[[i]]))
    data.frame(
      category = gene_tbl$category[[i]],
      gene = genes,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()

  main_gene_map <- main_gene_map %>%
    distinct(category, gene) %>%
    filter(gene %in% rownames(seu)) %>%
    mutate(category_clean = make_safe_name(category))

  main_genes <- unique(main_gene_map$gene)

  for (i in seq_len(nrow(sig_plot_tbl))) {
    score_col <- sig_plot_tbl$score_name[[i]]
    cat_name <- sig_plot_tbl$category[[i]]
    cat_safe <- sig_plot_tbl$category_clean[[i]]

    message_ts("Making signature violin plot for ", score_col)

    p_vln <- VlnPlot(
      object = seu,
      features = score_col,
      group.by = "seurat_clusters",
      pt.size = 0,
      combine = TRUE
    ) +
      retangio_theme_proj() +
      labs(
        title = paste0("Signature score — ", cat_name, " (broader)"),
        x = "Cluster",
        y = "Module score"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none"
      )

    retangio_save_plot(
      plot_obj = p_vln,
      filename_base = file.path(outdir, paste0("08i_signature_violin_", cat_safe)),
      width = 7,
      height = 6
    )
  }

  DefaultAssay(seu) <- "RNA"

  for (i in seq_len(nrow(sig_plot_tbl))) {
    score_col <- sig_plot_tbl$score_name[[i]]
    cat_name <- sig_plot_tbl$category[[i]]
    cat_safe <- sig_plot_tbl$category_clean[[i]]

    message_ts("Making signature feature plot for ", score_col)

    p_feat <- FeaturePlot(
      object = seu,
      features = score_col,
      reduction = reduction_to_use,
      order = TRUE,
      pt.size = 0.10,
      combine = TRUE
    ) +
      retangio_theme_proj() +
      labs(
        title = paste0("Signature score — ", cat_name, " (broader)"),
        x = "UMAP 1",
        y = "UMAP 2"
      )

    retangio_save_plot(
      plot_obj = p_feat,
      filename_base = file.path(outdir, paste0("08i_signature_feature_", cat_safe)),
      width = 7,
      height = 6
    )
  }

  heat_df <- cluster_scores %>%
    pivot_longer(
      cols = all_of(score_cols),
      names_to = "score_name",
      values_to = "mean_score"
    ) %>%
    left_join(sig_plot_tbl %>% select(score_name, category, category_clean), by = "score_name")

  heat_df$seurat_clusters <- factor(
    heat_df$seurat_clusters,
    levels = as.character(sort(as.numeric(unique(heat_df$seurat_clusters))))
  )

  heat_df$category <- factor(
    heat_df$category,
    levels = sig_plot_tbl$category
  )

  p_heat <- ggplot(heat_df, aes(x = category, y = seurat_clusters, fill = mean_score)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    retangio_theme_proj() +
    labs(
      title = "Cluster mean signature scores",
      x = "Signature category",
      y = "Cluster",
      fill = "Mean score"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )

  retangio_save_plot(
    plot_obj = p_heat,
    filename_base = file.path(outdir, "08i_signature_cluster_mean_heatmap"),
    width = 10,
    height = 6
  )

  DefaultAssay(seu) <- "RNA"

  for (g in main_genes) {
    message_ts("Making gene violin plot for ", g)

    p_gene_vln <- VlnPlot(
      object = seu,
      features = g,
      assay = "RNA",
      group.by = "seurat_clusters",
      pt.size = 0,
      combine = TRUE
    ) +
      retangio_theme_proj() +
      labs(
        title = paste0("Gene expression — ", g),
        x = "Cluster",
        y = "Expression level"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "none"
      )

    retangio_save_plot(
      plot_obj = p_gene_vln,
      filename_base = file.path(outdir, paste0("08i_gene_violin_", g)),
      width = 7,
      height = 6
    )
  }

  for (g in main_genes) {
    message_ts("Making gene feature plot for ", g)

    p_gene_feat <- FeaturePlot(
      object = seu,
      features = g,
      reduction = reduction_to_use,
      order = TRUE,
      pt.size = 0.10,
      combine = TRUE
    ) +
      retangio_theme_proj() +
      labs(
        title = paste0("Gene expression — ", g),
        x = "UMAP 1",
        y = "UMAP 2"
      )

    retangio_save_plot(
      plot_obj = p_gene_feat,
      filename_base = file.path(outdir, paste0("08i_gene_feature_", g)),
      width = 7,
      height = 6
    )
  }

  gene_dot_category_order_clean <- c(
    "arterial",
    "capillary_BRB",
    "venous",
    "tip",
    "proliferative"
  )

  gene_dot_category_labels <- c(
    "arterial" = "arterial",
    "capillary_BRB" = "capillary/BRB",
    "venous" = "venous",
    "tip" = "tip",
    "proliferative" = "proliferative"
  )

  gene_dot_map <- main_gene_map %>%
    mutate(category_clean = make_safe_name(category)) %>%
    filter(category_clean %in% gene_dot_category_order_clean) %>%
    mutate(
      category_clean = factor(category_clean, levels = gene_dot_category_order_clean),
      category_label = gene_dot_category_labels[as.character(category_clean)]
    ) %>%
    arrange(category_clean, gene)

  gene_order <- gene_dot_map$gene
  category_lookup_clean <- setNames(as.character(gene_dot_map$category_clean), gene_dot_map$gene)
  category_lookup_label <- setNames(as.character(gene_dot_map$category_label), gene_dot_map$gene)

  expr <- GetAssayData(seu, assay = "RNA", layer = "data")[gene_order, , drop = FALSE]
  meta <- seu@meta.data
  clusters <- levels(factor(
    as.character(seu$seurat_clusters),
    levels = as.character(sort(unique(as.numeric(as.character(seu$seurat_clusters)))))
  ))

  avg_exp <- sapply(clusters, function(cl) {
    cells <- rownames(meta)[as.character(meta$seurat_clusters) == cl]
    Matrix::rowMeans(expr[, cells, drop = FALSE])
  })

  pct_exp <- sapply(clusters, function(cl) {
    cells <- rownames(meta)[as.character(meta$seurat_clusters) == cl]
    Matrix::rowMeans(expr[, cells, drop = FALSE] > 0) * 100
  })

  avg_df <- as.data.frame(avg_exp)
  avg_df$gene <- rownames(avg_df)
  pct_df <- as.data.frame(pct_exp)
  pct_df$gene <- rownames(pct_df)

  plot_df <- avg_df %>%
    pivot_longer(-gene, names_to = "cluster", values_to = "avg_exp") %>%
    left_join(
      pct_df %>% pivot_longer(-gene, names_to = "cluster", values_to = "pct_exp"),
      by = c("gene", "cluster")
    ) %>%
    mutate(
      category_clean = category_lookup_clean[gene],
      category_label = category_lookup_label[gene],
      category_clean = factor(category_clean, levels = gene_dot_category_order_clean),
      category_label = factor(
        category_label,
        levels = gene_dot_category_labels[gene_dot_category_order_clean]
      ),
      gene = factor(gene, levels = unique(gene_order)),
      cluster = factor(cluster, levels = clusters)
    )

  p_dot <- ggplot(plot_df, aes(x = gene, y = cluster)) +
    geom_point(aes(size = pct_exp, color = avg_exp)) +
    facet_grid(. ~ category_label, scales = "free_x", space = "free_x") +
    scale_size(name = "Percent Expressed", range = c(0.5, 8)) +
    scale_color_gradient2(
      name = "Average Expression",
      low = "#D9D9D9",
      mid = "#B39DDB",
      high = "#1F00FF",
      midpoint = 0
    ) +
    retangio_theme_proj() +
    labs(
      title = "Main individual genes by cluster",
      x = "Genes",
      y = "Cluster"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.spacing.x = unit(0.6, "lines"),
      strip.background = element_rect(color = "black", fill = "white", linewidth = 1.2),
      strip.text = element_text(face = "bold"),
      legend.position = "right"
    )

  retangio_save_plot(
    plot_obj = p_dot,
    filename_base = file.path(outdir, "08i_gene_dotplot_main_genes"),
    width = 16,
    height = 8
  )

  message_ts("Stage 08i_plot_signatures_and_genes_ec_after_soupx_unintegrated complete")
}

main()
