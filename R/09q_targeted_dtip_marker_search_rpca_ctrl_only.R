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
  markers
}

main <- function() {
  log_message("Starting stage 09q_targeted_dtip_marker_search_rpca_ctrl_only")

  future::plan("sequential")
  options(future.globals.maxSize = 8 * 1024^3)

  input_rds <- "results/objects/09o_rpca_tipcell_ctrl_only_resolution_sweep.rds"
  out_dir_markers <- "results/markers"
  out_dir_review <- "results/review"

  dir.create(out_dir_markers, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_dir_review, recursive = TRUE, showWarnings = FALSE)

  seu <- readRDS(input_rds)
  log_message("Loaded object: ", input_rds)

  cluster_col <- "tip_ctrl_res_0_20"
  if (!cluster_col %in% colnames(seu@meta.data)) {
    stop("Missing cluster column: ", cluster_col)
  }

  clusters_keep <- c("0", "1", "2")

  seu <- subset(
    seu,
    cells = colnames(seu)[as.character(seu@meta.data[[cluster_col]]) %in% clusters_keep]
  )
  log_message("Cells retained: ", ncol(seu))

  seu$tip_ctrl_cluster_targeted <- factor(
    as.character(seu@meta.data[[cluster_col]]),
    levels = clusters_keep
  )

  seu$comparison_0_vs_1 <- ifelse(
    as.character(seu$tip_ctrl_cluster_targeted) == "0",
    "cluster0",
    ifelse(as.character(seu$tip_ctrl_cluster_targeted) == "1", "cluster1", NA)
  )

  seu$comparison_0_vs_12 <- ifelse(
    as.character(seu$tip_ctrl_cluster_targeted) == "0",
    "cluster0",
    "clusters1_2"
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

  seu_01 <- subset(seu, subset = !is.na(comparison_0_vs_1))
  log_message("Cells in cluster0 vs cluster1 comparison: ", ncol(seu_01))

  m1 <- run_one_comparison(
    seu = seu_01,
    group_col = "comparison_0_vs_1",
    ident_1 = "cluster0",
    ident_2 = "cluster1",
    label_1 = "cluster0",
    label_2 = "cluster1",
    out_csv = file.path(out_dir_markers, "09q_cluster0_vs_cluster1_markers.csv")
  )
  log_message("Saved cluster 0 vs cluster 1 comparison")

  m2 <- run_one_comparison(
    seu = seu,
    group_col = "comparison_0_vs_12",
    ident_1 = "cluster0",
    ident_2 = "clusters1_2",
    label_1 = "cluster0",
    label_2 = "clusters1_2",
    out_csv = file.path(out_dir_markers, "09q_cluster0_vs_clusters1_2_markers.csv")
  )
  log_message("Saved cluster 0 vs clusters 1+2 comparison")

  topN <- function(df, n = 30) {
    df %>%
      filter(!is.na(avg_log2FC), !is.na(p_val_adj)) %>%
      slice_max(order_by = avg_log2FC, n = n, with_ties = FALSE)
  }

  top30_m1 <- topN(m1, 30)
  top30_m2 <- topN(m2, 30)

  write_csv(
    top30_m1,
    file.path(out_dir_review, "09q_cluster0_vs_cluster1_top30_markers.csv")
  )
  write_csv(
    top30_m2,
    file.path(out_dir_review, "09q_cluster0_vs_clusters1_2_top30_markers.csv")
  )

  summary_tbl <- bind_rows(
    m1 %>% slice_max(order_by = avg_log2FC, n = 15, with_ties = FALSE),
    m2 %>% slice_max(order_by = avg_log2FC, n = 15, with_ties = FALSE)
  ) %>%
    select(comparison, gene, avg_log2FC, pct.1, pct.2, p_val_adj)

  write_csv(
    summary_tbl,
    file.path(out_dir_review, "09q_ctrl_only_dtip_candidate_marker_summary_top15.csv")
  )

  log_message("Stage 09q_targeted_dtip_marker_search_rpca_ctrl_only complete")
}

main()
