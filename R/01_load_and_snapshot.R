source("R/00_utils.R")

main <- function() {
  cfg <- read_config()
  outdir <- cfg$project$output_dir
  ensure_dir(outdir)
  set_seed_from_config(cfg)

  in_obj <- cfg$input$seurat_rds
  if (!file.exists(in_obj)) {
    fail_fast("Input Seurat object not found: {in_obj}")
  }

  log_info("Loading Seurat object from {in_obj}")
  seu <- load_object(in_obj)

  save_object(seu, file.path(outdir, "01_loaded_raw"), cfg, object_type = "seurat")

  meta_cols <- c(
    "orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo",
    "doublettype", "doubletscore", cfg$integration$final_cluster_column
  )
  present <- intersect(meta_cols, colnames(seu@meta.data))
  readr::write_csv(seu@meta.data[, present, drop = FALSE], file.path(outdir, "01_metadata_snapshot.csv"))

  save_config_snapshot(cfg, outdir)
  save_session_info(outdir)
  save_package_versions(outdir)
  run_io_sanity_checks(cfg, outdir)

  log_info("01_load_and_snapshot complete")
}

if (sys.nframe() == 0) {
  main()
}
