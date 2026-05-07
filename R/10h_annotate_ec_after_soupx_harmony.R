#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(tibble)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

main <- function() {
  input_rds <- "results/objects/10f_ec_after_soupx_harmony_signature_scored.rds"

  output_annotation_table <- "results/annotation/10h_ec_after_soupx_harmony_cluster_annotation_table.csv"
  output_metadata_table <- "results/annotation/10h_ec_after_soupx_harmony_annotated_metadata.csv"
  output_rds <- "results/objects/10h_ec_after_soupx_harmony_annotated.rds"

  dir.create("results/annotation", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 10h_annotate_ec_after_soupx_harmony")

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    stop("seurat_clusters not found in metadata.")
  }

  annotation_tbl <- tibble(
    seurat_clusters = c("0","1","2","3","4","5","6","7","8","9"),
    proposed_label = c(
      "capillary_BRB_venous_biased",
      "tip_angiogenic",
      "arterial_side_capillary_BRB_homeostatic",
      "capillary_BRB",
      "proliferative_S_phase",
      "atypical_capillary_BRB_stressed_mixed",
      "specialized_Aqp1_Six3_like_EC",
      "proliferative_G2M",
      "inflammatory_IFN_EC",
      "rare_specialized_EC"
    ),
    my_label = c(
      "Venous x Capillary",
      "Tip cell",
      "Arterial x Capillary",
      "Capillary",
      "Proliferative 1",
      "Uncertain",
      "EC 1",
      "Proliferative 2",
      "Inflammatory",
      "EC 2"
    )
  )

  observed_clusters <- sort(unique(as.character(seu$seurat_clusters)))
  expected_clusters <- sort(annotation_tbl$seurat_clusters)

  if (!identical(observed_clusters, expected_clusters)) {
    stop(
      paste0(
        "Cluster mismatch between object and annotation table.\n",
        "Observed: ", paste(observed_clusters, collapse = ", "), "\n",
        "Expected: ", paste(expected_clusters, collapse = ", ")
      )
    )
  }

  write_csv(annotation_tbl, output_annotation_table)
  message_ts("Saved annotation table:", output_annotation_table)

  map_proposed <- setNames(annotation_tbl$proposed_label, annotation_tbl$seurat_clusters)
  map_my <- setNames(annotation_tbl$my_label, annotation_tbl$seurat_clusters)

  seu$cluster_id <- as.character(seu$seurat_clusters)
  seu$annotation_proposed <- unname(map_proposed[seu$cluster_id])
  seu$annotation_my <- unname(map_my[seu$cluster_id])

  cluster_order_num <- sort(unique(as.numeric(annotation_tbl$seurat_clusters)))
  cluster_order_chr <- as.character(cluster_order_num)

  my_label_levels <- annotation_tbl$my_label[match(cluster_order_chr, annotation_tbl$seurat_clusters)]
  proposed_levels <- annotation_tbl$proposed_label[match(cluster_order_chr, annotation_tbl$seurat_clusters)]

  seu$annotation_my <- factor(seu$annotation_my, levels = my_label_levels)
  seu$annotation_proposed <- factor(seu$annotation_proposed, levels = proposed_levels)

  Idents(seu) <- "annotation_my"

  annotated_metadata <- seu@meta.data %>%
    rownames_to_column("cell_id")

  write_csv(annotated_metadata, output_metadata_table)
  message_ts("Saved annotated metadata:", output_metadata_table)

  if (is.null(seu@misc)) {
    seu@misc <- list()
  }

  seu@misc$stage10h_annotation <- list(
    source_object = input_rds,
    annotation_table = output_annotation_table,
    active_ident = "annotation_my",
    run_datetime = as.character(Sys.time())
  )

  saveRDS(seu, output_rds)
  message_ts("Saved annotated object:", output_rds)
  message_ts("Number of annotated cells:", ncol(seu))
  message_ts("Annotation labels:", paste(levels(seu$annotation_my), collapse = ", "))
  message_ts("Stage 10h_annotate_ec_after_soupx_harmony complete")
}

main()
