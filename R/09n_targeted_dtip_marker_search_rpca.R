#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(future)
})

log_message <- function(...) {
  msg <- paste0(...)
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", msg, "\n", sep = "")
}

run_one_comparison <- function(
  seu,
  group_col,
  ident_1,
  ident_2,
  label_1,
  label_2,
  out_csv
) {
  Idents(seu) <- seu@meta.data[[group_col]]

  markers <- FindMarkers(
    object = seu,
    ident.1 = ident_1,
    ident.2 = ident_2,
    assay = "RNA",
    test.use = "wilcox",
    only.pos = TRUE,
    min.pct = 0.10,
    logfc.threshold = 0.10,
    verbose = FALSE
  )

  markers <- markers %>%
    tibble::rownames_to_column("gene") %>%
    arrange(desc(avg_log2FC), p_val_adj, desc(pct.1)) %>%
    mutate(
      comparison = paste0(label_1, "_vs_", label_2),
      group_1 = label_1,
      group_2 = label_2
    ) %>%
    select(comparison, group_1, group_2, gene, everything())

  write_csv(markers, out_csv)

  return(markers)
}

main <- function() {
  log_message("Starting stage 09n_targeted_dtip_marker_search_rpca")

  future::plan("sequential")
  options(future.globals.maxSize = 8 * 1024^3)

  input_rds <- "results/objects/09k_rpca_tipcell_resolution_sweep.rds"
  out_dir_markers <- "results/markers"
  out_dir_review <- "results/review"

  dir.create(out_dir_markers, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir_review, recursive = TRUE, showWarnings = FALSE)

  seu <- readRDS(input_rds)
  log_message("Loaded object: ", input_rds)

  cluster_col <- "tip_res_0_20"
  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop("Missing cluster column: ", cluster_col)
  }

  clusters_keep <- c("0", "1", "2", "3")
  clusters_exclude <- c("4", "5")

  seu <- subset(
    seu,
    cells = colnames(seu)[as.character(seu@meta.data[[cluster_col]]) %in% clusters_keep]
  )
  log_message("Cells retained after excluding clusters 4 and 5: ", ncol(seu))

  seu$tip_cluster_targeted <- factor(
    as.character(seu@meta.data[[cluster_col]]),
    levels = clusters_keep
  )

  seu$comparison_1_vs_rest <- ifelse(
    as.character(seu$tip_cluster_targeted) == "1",
    "cluster1",
    "other_0_2_3"
  )

  seu$comparison_2_vs_rest <- ifelse(
    as.character(seu$tip_cluster_targeted) == "2",
    "cluster2",
    "other_0_1_3"
  )

  seu$comparison_12_vs_03 <- ifelse(
    as.character(seu$tip_cluster_targeted) %in% c("1", "2"),
    "clusters1_2",
    "clusters0_3"
  )

  DefaultAssay(seu) <- "RNA"
  log_message("Joining RNA layers")
  seu <- JoinLayers(seu, assay = "RNA")

  log_message("Normalizing joined RNA assay")
  seu <- NormalizeData(
    object = seu,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  m1 <- run_one_comparison(
    seu = seu,
    group_col = "comparison_1_vs_rest",
    ident_1 = "cluster1",
    ident_2 = "other_0_2_3",
    label_1 = "cluster1",
    label_2 = "other_0_2_3",
    out_csv = file.path(out_dir_markers, "09n_cluster1_vs_0_2_3_markers.csv")
  )
  log_message("Saved cluster 1 comparison")

  m2 <- run_one_comparison(
    seu = seu,
    group_col = "comparison_2_vs_rest",
    ident_1 = "cluster2",
    ident_2 = "other_0_1_3",
    label_1 = "cluster2",
    label_2 = "other_0_1_3",
    out_csv = file.path(out_dir_markers, "09n_cluster2_vs_0_1_3_markers.csv")
  )
  log_message("Saved cluster 2 comparison")

  m3 <- run_one_comparison(
    seu = seu,
    group_col = "comparison_12_vs_03",
    ident_1 = "clusters1_2",
    ident_2 = "clusters0_3",
    label_1 = "clusters1_2",
    label_2 = "clusters0_3",
    out_csv = file.path(out_dir_markers, "09n_clusters1_2_vs_0_3_markers.csv")
  )
  log_message("Saved combined comparison")

  topN <- function(df, n = 30) {
    df %>%
      filter(!is.na(avg_log2FC), !is.na(p_val_adj)) %>%
      slice_max(order_by = avg_log2FC, n = n, with_ties = FALSE)
  }

  write_csv(
    topN(m1, 30),
    file.path(out_dir_review, "09n_cluster1_vs_0_2_3_top30_markers.csv")
  )
  write_csv(
    topN(m2, 30),
    file.path(out_dir_review, "09n_cluster2_vs_0_1_3_top30_markers.csv")
  )
  write_csv(
    topN(m3, 30),
    file.path(out_dir_review, "09n_clusters1_2_vs_0_3_top30_markers.csv")
  )

  summary_tbl <- bind_rows(
    m1 %>% slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE),
    m2 %>% slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE),
    m3 %>% slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE)
  ) %>%
    select(comparison, gene, avg_log2FC, pct.1, pct.2, p_val_adj)

  write_csv(
    summary_tbl,
    file.path(out_dir_review, "09n_dtip_candidate_marker_summary_top10.csv")
  )

  log_message("Stage 09n_targeted_dtip_marker_search_rpca complete")
}

main()
