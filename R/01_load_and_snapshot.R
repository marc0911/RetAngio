source("R/00_utils.R")
suppressPackageStartupMessages({
  library(qs)
})

main <- function() {
  cfg <- read_config()
  outdir <- cfg$project$output_dir
  ensure_dir(outdir)
  set_seed_from_config(cfg)

  in_rds <- cfg$input$seurat_rds
  if (!file.exists(in_rds)) {
    fail_fast("Input Seurat object not found: {in_rds}")
  }

  log_info("Loading Seurat object from {in_rds}")
  seu <- readRDS(in_rds)

  saveRDS(seu, file.path(outdir, "01_loaded_raw.rds"))
  qs::qsave(seu, file.path(outdir, "01_loaded_raw.qs"), preset = "high")

  meta_cols <- c(
    "orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo",
    "doublettype", "doubletscore", cfg$integration$final_cluster_column
  )
  present <- intersect(meta_cols, colnames(seu@meta.data))
  readr::write_csv(seu@meta.data[, present, drop = FALSE], file.path(outdir, "01_metadata_snapshot.csv"))

  save_config_snapshot(cfg, outdir)
  save_session_info(outdir)
  save_package_versions(outdir)

  log_info("01_load_and_snapshot complete")
}

if (sys.nframe() == 0) {
  main()
}
