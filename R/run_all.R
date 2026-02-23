source("R/00_utils.R")

main <- function() {
  scripts <- c(
    "R/01_load_and_snapshot.R",
    "R/02_remove_contaminants_and_reintegrate.R",
    "R/03_build_BASE_object.R",
    "R/04_QC_and_contamination_checks.R",
    "R/05_annotation_program_scoring.R",
    "R/06_DGE_cluster_markers_and_pseudobulk.R"
  )

  for (s in scripts) {
    log_info("Running {s}")
    e <- new.env(parent = globalenv())
    sys.source(s, envir = e)
    e$main()
  }

  cfg <- read_config()
  outdir <- cfg$project$output_dir
  seu <- load_object(file.path(outdir, "05_annotated"), cfg)

  check_contaminants_removed(seu, "final_cluster", cfg$integration$contaminants)
  check_required_layers(seu, "RNA")
  assert_ident_matches(seu, "final_cluster")

  qc_path <- file.path(outdir, "06_pseudobulk_qc_table.csv")
  if (!file.exists(qc_path)) fail_fast("Missing pseudobulk QC table: {qc_path}")
  qc_tbl <- readr::read_csv(qc_path, show_col_types = FALSE)
  if (!(any(qc_tbl$status == "tested") || all(qc_tbl$status == "skipped"))) {
    fail_fast("Pseudobulk QC invalid: must have tested clusters or explicitly skipped all")
  }

  log_info("run_all complete")
}

if (sys.nframe() == 0) {
  main()
}
