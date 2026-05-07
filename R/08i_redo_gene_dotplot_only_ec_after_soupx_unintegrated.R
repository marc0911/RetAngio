#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readxl)
  library(ggplot2)
  library(tidyr)
  library(Matrix)
})

source("R/plotting_helpers.R")

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

clean_gene_vector <- function(x) {
  x <- unlist(strsplit(x, ","))
  x <- trimws(x)
  x <- x[nchar(x) > 0]
  unique(x)
}

make_safe_name <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

main <- function() {
  input_rds <- "results/objects/08h_ec_after_soupx_unintegrated_signature_scored.rds"
  input_xlsx <- "reference/RetAngio_annotation_master_v2.xlsx"
  outdir <- "plots/08i_ec_after_soupx_unintegrated_signatures_and_genes"

  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting 08i redo gene dotplot only")

  seu <- readRDS(input_rds)
  sig_tbl <- readxl::read_excel(input_xlsx, sheet = "category_signatures")

  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    seu$seurat_clusters <- as.character(Idents(seu))
  }
  seu$seurat_clusters <- factor(
    as.character(seu$seurat_clusters),
    levels = as.character(sort(unique(as.numeric(as.character(seu$seurat_clusters)))))
  )

  if (!"RNA" %in% Assays(seu)) {
    stop("RNA assay not found in object.")
  }

  DefaultAssay(seu) <- "RNA"
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  seu <- NormalizeData(seu, verbose = FALSE)

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

  gene_tbl <- sig_tbl %>%
    transmute(
      category_raw = category,
      category_clean = make_safe_name(category),
      minimal = recommended_minimal_signature
    ) %>%
    filter(category_clean %in% gene_dot_category_order_clean)

  main_gene_map <- lapply(seq_len(nrow(gene_tbl)), function(i) {
    genes <- clean_gene_vector(as.character(gene_tbl$minimal[[i]]))
    data.frame(
      category_raw = gene_tbl$category_raw[[i]],
      category_clean = gene_tbl$category_clean[[i]],
      gene = genes,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()

  main_gene_map <- main_gene_map %>%
    distinct(category_raw, category_clean, gene) %>%
    filter(gene %in% rownames(seu)) %>%
    mutate(
      category_clean = factor(category_clean, levels = gene_dot_category_order_clean),
      category_label = gene_dot_category_labels[as.character(category_clean)]
    ) %>%
    arrange(category_clean, gene)

  gene_order <- main_gene_map$gene
  category_lookup_clean <- setNames(as.character(main_gene_map$category_clean), main_gene_map$gene)
  category_lookup_label <- setNames(as.character(main_gene_map$category_label), main_gene_map$gene)

  expr <- GetAssayData(seu, assay = "RNA", layer = "data")[gene_order, , drop = FALSE]
  meta <- seu@meta.data
  clusters <- levels(meta$seurat_clusters)

  avg_exp <- sapply(clusters, function(cl) {
    cells <- rownames(meta)[meta$seurat_clusters == cl]
    Matrix::rowMeans(expr[, cells, drop = FALSE])
  })

  pct_exp <- sapply(clusters, function(cl) {
    cells <- rownames(meta)[meta$seurat_clusters == cl]
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

  message_ts("Rewrote dotplot: ", file.path(outdir, "08i_gene_dotplot_main_genes"))
  message_ts("08i redo gene dotplot only complete")
}

main()
