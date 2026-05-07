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
  input_rds <- "results/objects/09f_ec_after_soupx_rpca_signature_scored.rds"

  output_annotation_table <- "results/annotation/09h_ec_after_soupx_rpca_cluster_annotation_table.csv"
  output_metadata_table <- "results/annotation/09h_ec_after_soupx_rpca_annotated_metadata.csv"
  output_rds <- "results/objects/09h_ec_after_soupx_rpca_annotated.rds"

  dir.create("results/annotation", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 09h_annotate_ec_after_soupx_rpca")

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    stop("seurat_clusters not found in metadata.")
  }

  annotation_tbl <- tibble(
    seurat_clusters = c("0","1","2","3","4","5","6","7","8","9","10","11","12"),
    proposed_label = c(
      "venous_side_capillary_BRB",
      "capillary_BRB",
      "capillary_BRB_specialized",
      "arterial_side_capillary_BRB_Cxcl12_high",
      "proliferative_S_phase",
      "tip_angiogenic",
      "arterial",
      "specialized_Aqp1_Six3_like_EC",
      "proliferative_G2M",
      "atypical_capillary_BRB_stressed_mixed",
      "proliferative_mitotic",
      "inflammatory_IFN_EC",
      "rare_specialized_EC"
    ),
    my_label = c(
      "Venous x capillary",
      "Capillary 1",
      "Capillary 2",
      "Arterial x capillary",
      "Proliferative 1",
      "Tip cell",
      "Arterial",
      "EC1",
      "Proliferative 2",
      "Uncertain",
      "Proliferative 3",
      "Inflammatory EC",
      "EC2"
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

  my_label_levels <- annotation_tbl$my_label[match(sort(unique(as.numeric(annotation_tbl$seurat_clusters))), as.numeric(annotation_tbl$seurat_clusters))]
  seu$annotation_my <- factor(seu$annotation_my, levels = my_label_levels)

  proposed_levels <- annotation_tbl$proposed_label[match(sort(unique(as.numeric(annotation_tbl$seurat_clusters))), as.numeric(annotation_tbl$seurat_clusters))]
  seu$annotation_proposed <- factor(seu$annotation_proposed, levels = proposed_levels)

  Idents(seu) <- "annotation_my"

  annotated_metadata <- seu@meta.data %>%
    rownames_to_column("cell_id")

  write_csv(annotated_metadata, output_metadata_table)
  message_ts("Saved annotated metadata:", output_metadata_table)

  if (is.null(seu@misc)) {
    seu@misc <- list()
  }

  seu@misc$stage09h_annotation <- list(
    source_object = input_rds,
    annotation_table = output_annotation_table,
    active_ident = "annotation_my",
    run_datetime = as.character(Sys.time())
  )

  saveRDS(seu, output_rds)
  message_ts("Saved annotated object:", output_rds)
  message_ts("Number of annotated cells:", ncol(seu))
  message_ts("Annotation labels:", paste(levels(seu$annotation_my), collapse = ", "))
  message_ts("Stage 09h_annotate_ec_after_soupx_rpca complete")
}

main()
