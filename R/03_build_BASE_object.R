source("R/00_utils.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(lobstr)
  library(readr)
})

main <- function() {
  cfg <- read_config()
  outdir <- cfg$project$output_dir
  ensure_dir(outdir)

  seu <- readRDS(file.path(outdir, "02_clean_reintegrated.rds"))
  if (!"final_cluster" %in% colnames(seu@meta.data)) {
    fail_fast("final_cluster missing; run script 02 first")
  }
  Idents(seu) <- "final_cluster"

  if (!"RNA" %in% names(seu@assays)) fail_fast("RNA assay missing")
  if (!"SCT" %in% names(seu@assays)) fail_fast("SCT assay missing")

  # Ensure RNA has counts + data layers for marker and DEG use.
  if (!"data" %in% SeuratObject::Layers(seu[["RNA"]])) {
    seu <- NormalizeData(seu, assay = "RNA", verbose = FALSE)
  }
  check_required_layers(seu, "RNA")

  # Keep counts only in SCT for memory-safe BASE snapshot.
  sct_counts <- get_layer_data(seu, assay = "SCT", layer = "counts", fallback_slot = "counts")
  seu[["SCT"]] <- CreateAssayObject(counts = sct_counts)

  essential <- c(
    cfg$input$sample_column,
    cfg$input$condition_column,
    cfg$input$timepoint_column,
    "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "final_cluster"
  )
  keep_cols <- intersect(essential, colnames(seu@meta.data))
  seu@meta.data <- seu@meta.data[, keep_cols, drop = FALSE]

  assert_ident_matches(seu, "final_cluster")
  assert_reductions_present(seu)
  check_contaminants_removed(seu, "final_cluster", cfg$integration$contaminants)

  out_rds <- file.path(outdir, "03_BASE_object.rds")
  saveRDS(seu, out_rds, compress = "xz")

  meta <- tibble::tibble(
    object = "03_BASE_object.rds",
    object_size_bytes = lobstr::obj_size(seu),
    cells = ncol(seu),
    genes = nrow(seu)
  )
  write_csv(meta, file.path(outdir, "03_BASE_object_summary.csv"))
  write_checksum(out_rds, file.path(outdir, "03_BASE_object_checksum.csv"))

  save_session_info(outdir)
  save_config_snapshot(cfg, outdir)
  log_info("03_build_BASE_object complete")
}

if (sys.nframe() == 0) {
  main()
}
