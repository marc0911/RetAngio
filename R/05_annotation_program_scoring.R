source("R/00_utils.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readxl)
  library(tidyr)
  library(readr)
})

read_zarkada_sets <- function(path, sheet) {
  if (!file.exists(path)) fail_fast("Missing Zarkada Excel file: {path}")
  raw <- readxl::read_excel(path, sheet = sheet)
  req <- c("program", "gene")
  if (!all(req %in% colnames(raw))) {
    fail_fast("Zarkada sheet must include columns: program, gene")
  }
  split(raw$gene, raw$program)
}

merge_gene_sets <- function(gene_sets, merge_cfg) {
  if (is.null(merge_cfg)) return(gene_sets)
  for (new_name in names(merge_cfg)) {
    parts <- merge_cfg[[new_name]]
    genes <- unique(unlist(gene_sets[parts]))
    gene_sets[[new_name]] <- genes
  }
  gene_sets
}

main <- function() {
  cfg <- read_config()
  outdir <- cfg$project$output_dir
  ensure_dir(outdir)

  seu <- load_object(file.path(outdir, "03_BASE_object"), cfg)
  gene_sets <- read_zarkada_sets(cfg$input$zarkada_excel, cfg$input$zarkada_sheet)
  gene_sets <- merge_gene_sets(gene_sets, cfg$annotation$merge_sets)

  score_method <- tolower(cfg$annotation$score_method)
  if (score_method == "ucell") {
    # UCell robust rank-based scoring reference: PMID:34285779
    if (!requireNamespace("UCell", quietly = TRUE)) {
      fail_fast("UCell requested but package not installed")
    }
    seu <- UCell::AddModuleScore_UCell(seu, features = gene_sets, assay = "RNA")
    score_cols <- grep("_UCell$", colnames(seu@meta.data), value = TRUE)
  } else if (score_method == "addmodulescore") {
    seu <- AddModuleScore(seu, features = unname(gene_sets), assay = "RNA", name = names(gene_sets), seed = cfg$project$seed)
    score_cols <- grep(paste0("^(", paste(names(gene_sets), collapse = "|"), ")"), colnames(seu@meta.data), value = TRUE)
  } else {
    fail_fast("annotation.score_method must be 'ucell' or 'addmodulescore'")
  }

  Idents(seu) <- "final_cluster"

  summary_tbl <- seu@meta.data |>
    tibble::rownames_to_column("cell") |>
    group_by(final_cluster) |>
    summarise(across(all_of(score_cols), mean, na.rm = TRUE), .groups = "drop")
  write_csv(summary_tbl, file.path(outdir, "05_annotation_summary.csv"))

  save_object(seu, file.path(outdir, "05_annotated"), cfg, object_type = "seurat")
  log_info("05_annotation_program_scoring complete")
}

if (sys.nframe() == 0) {
  main()
}
