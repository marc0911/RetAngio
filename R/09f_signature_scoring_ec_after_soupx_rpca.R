#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(readxl)
})

message_ts <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "|", ..., "\n")
}

clean_gene_vector <- function(x) {
  x <- unlist(strsplit(x, ","))
  x <- trimws(x)
  x <- x[nchar(x) > 0]
  unique(x)
}

main <- function() {
  input_rds <- "results/objects/09c_ec_after_soupx_rpca_final.rds"
  input_xlsx <- "reference/RetAngio_annotation_master_v2.xlsx"

  output_rds <- "results/objects/09f_ec_after_soupx_rpca_signature_scored.rds"
  output_sig_used <- "results/signatures/09f_ec_after_soupx_rpca_signatures_used.csv"
  output_cluster_scores <- "results/signatures/09f_ec_after_soupx_rpca_cluster_mean_scores.csv"
  output_cell_scores <- "results/signatures/09f_ec_after_soupx_rpca_cell_scores.csv"

  dir.create("results/objects", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/signatures", recursive = TRUE, showWarnings = FALSE)

  message_ts("Starting stage 09f_signature_scoring_ec_after_soupx_rpca")

  seu <- readRDS(input_rds)
  message_ts("Loaded object:", input_rds)

  if (!file.exists(input_xlsx)) {
    stop(paste("Workbook not found:", input_xlsx))
  }

  sig_tbl <- readxl::read_excel(input_xlsx, sheet = "category_signatures")
  message_ts("Loaded workbook sheet: category_signatures")

  if (!all(c("category", "recommended_minimal_signature", "recommended_broader_signature") %in% colnames(sig_tbl))) {
    stop("Workbook sheet category_signatures is missing required columns.")
  }

  if (!"RNA" %in% Assays(seu)) {
    stop("RNA assay not found in object.")
  }

  DefaultAssay(seu) <- "RNA"
  seu[["RNA"]] <- JoinLayers(seu[["RNA"]])
  seu <- NormalizeData(seu, verbose = FALSE)

  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    seu$seurat_clusters <- as.character(Idents(seu))
  }

  genes_present <- rownames(seu)

  sig_rows <- list()
  sig_id_counter <- 1

  for (i in seq_len(nrow(sig_tbl))) {
    cat_name <- as.character(sig_tbl$category[[i]])

    min_genes <- clean_gene_vector(as.character(sig_tbl$recommended_minimal_signature[[i]]))
    broad_genes <- clean_gene_vector(as.character(sig_tbl$recommended_broader_signature[[i]]))

    min_present <- min_genes[min_genes %in% genes_present]
    broad_present <- broad_genes[broad_genes %in% genes_present]

    sig_rows[[length(sig_rows) + 1]] <- data.frame(
      category = cat_name,
      signature_set = "minimal",
      score_name = paste0("sig_", sig_id_counter, "_", gsub("[^A-Za-z0-9]+", "_", cat_name), "_minimal"),
      n_input_genes = length(min_genes),
      n_present_genes = length(min_present),
      genes_used = paste(min_present, collapse = ", "),
      stringsAsFactors = FALSE
    )
    sig_id_counter <- sig_id_counter + 1

    sig_rows[[length(sig_rows) + 1]] <- data.frame(
      category = cat_name,
      signature_set = "broader",
      score_name = paste0("sig_", sig_id_counter, "_", gsub("[^A-Za-z0-9]+", "_", cat_name), "_broader"),
      n_input_genes = length(broad_genes),
      n_present_genes = length(broad_present),
      genes_used = paste(broad_present, collapse = ", "),
      stringsAsFactors = FALSE
    )
    sig_id_counter <- sig_id_counter + 1
  }

  sig_used <- bind_rows(sig_rows)

  if (any(sig_used$n_present_genes == 0)) {
    bad <- sig_used %>% filter(n_present_genes == 0)
    stop(
      paste0(
        "Some signatures have zero genes present in the object: ",
        paste(bad$score_name, collapse = ", ")
      )
    )
  }

  write_csv(sig_used, output_sig_used)
  message_ts("Saved signature definitions used:", output_sig_used)

  for (i in seq_len(nrow(sig_used))) {
    score_name <- sig_used$score_name[[i]]
    genes_used <- clean_gene_vector(sig_used$genes_used[[i]])

    message_ts("Scoring signature: ", score_name, " (", length(genes_used), " genes)")

    seu <- AddModuleScore(
      object = seu,
      features = list(genes_used),
      assay = "RNA",
      name = paste0(score_name, "_"),
      search = FALSE,
      seed = 1
    )

    generated_col <- paste0(score_name, "_1")
    if (!generated_col %in% colnames(seu@meta.data)) {
      stop(paste("Expected score column not found after AddModuleScore:", generated_col))
    }

    seu@meta.data[[score_name]] <- seu@meta.data[[generated_col]]
    seu@meta.data[[generated_col]] <- NULL
  }

  score_cols <- sig_used$score_name

  cell_scores <- seu@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    select(cell_id, seurat_clusters, orig.ident, any_of(score_cols))

  write_csv(cell_scores, output_cell_scores)
  message_ts("Saved cell-level signature scores:", output_cell_scores)

  cluster_scores <- cell_scores %>%
    group_by(seurat_clusters) %>%
    summarise(across(all_of(score_cols), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
    arrange(as.numeric(as.character(seurat_clusters)))

  write_csv(cluster_scores, output_cluster_scores)
  message_ts("Saved cluster mean signature scores:", output_cluster_scores)

  if (is.null(seu@misc)) {
    seu@misc <- list()
  }

  seu@misc$stage09f_signature_scoring <- list(
    source_object = input_rds,
    signature_workbook = input_xlsx,
    signature_sheet = "category_signatures",
    score_columns = score_cols,
    run_datetime = as.character(Sys.time())
  )

  saveRDS(seu, output_rds)
  message_ts("Saved signature-scored object:", output_rds)
  message_ts("Stage 09f_signature_scoring_ec_after_soupx_rpca complete")
}

main()
