#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_rds <- "results/objects/10c_ec_after_soupx_harmony_final.rds"
  input_markers <- "results/markers/10d_ec_after_soupx_harmony_findallmarkers.csv"

  out_cluster_sizes <- "results/review/10e_ec_after_soupx_harmony_cluster_sizes.csv"
  out_top20 <- "results/review/10e_ec_after_soupx_harmony_top20_markers.csv"
  out_panel <- "results/review/10e_ec_after_soupx_harmony_marker_panel_avgexp.csv"

  dir.create("results/review", recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 10e_review_tables_ec_after_soupx_harmony")

  seu <- readRDS(input_rds)
  markers <- read_csv(input_markers, show_col_types = FALSE)

  message_ts("Loaded object:", input_rds)
  message_ts("Loaded markers:", input_markers)

  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    seu$seurat_clusters <- as.character(Idents(seu))
  }

  cluster_sizes <- seu@meta.data %>%
    count(seurat_clusters, name = "n_cells") %>%
    mutate(freq = n_cells / sum(n_cells)) %>%
    arrange(as.numeric(as.character(seurat_clusters)))

  write_csv(cluster_sizes, out_cluster_sizes)

  top20 <- markers %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = 20) %>%
    ungroup()

  write_csv(top20, out_top20)

  DefaultAssay(seu) <- "RNA"
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  seu <- NormalizeData(seu, verbose = FALSE)

  panel_genes <- c(
    "Sox17", "Nrp2", "Des", "Sag", "Glul", "Rgs5", "Cd34", "Snap25",
    "Gja5", "Vcam1", "Emcn", "Tek", "Elavl4", "Gja4", "C1qa", "Kit",
    "Kdr", "Pde6b", "Col1a2", "Rho", "Vwf", "Slco1a4", "Tyrobp", "Nr2f2",
    "Rlbp1", "Hbb-bs", "Adm", "Mki67", "Efnb2", "Angpt2", "Klf2", "Cspg4",
    "Pclaf", "Gnat1", "Dcn", "Lum", "Syt1", "Ptprb", "Lyz2", "Hba-a1",
    "Hba-a2", "Col1a1", "Top2a", "Pecam1", "Birc5", "Rbfox3", "Esm1",
    "Slc1a3", "Cldn5", "Notch3", "Aif1", "Mog", "Prph2", "Aqp4", "Pdgfrb",
    "Mbp", "Apln", "Plp1", "Bmx"
  )

  panel_genes <- panel_genes[panel_genes %in% rownames(seu)]

  avg_exp <- AverageExpression(
    object = seu,
    assays = "RNA",
    features = panel_genes,
    group.by = "seurat_clusters",
    slot = "data",
    verbose = FALSE
  )$RNA

  avg_exp_df <- as.data.frame(avg_exp)
  avg_exp_df$gene <- rownames(avg_exp_df)
  avg_exp_df <- avg_exp_df %>% relocate(gene)

  write_csv(avg_exp_df, out_panel)

  message_ts("Saved cluster sizes:", out_cluster_sizes)
  message_ts("Saved top20 markers:", out_top20)
  message_ts("Saved marker panel average expression:", out_panel)
  message_ts("Stage 10e_review_tables_ec_after_soupx_harmony complete")
}

main()
