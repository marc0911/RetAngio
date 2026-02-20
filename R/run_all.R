source("R/00_utils.R")

args <- parse_args()
cfg <- load_config(args$config)
set.seed(cfg$seed)
init_dirs(cfg)

scripts <- c(
  "R/01_load_and_snapshot.R",
  "R/02_remove_contaminants_and_recluster.R",
  "R/03_build_BASE_object.R",
  "R/04_QC_and_contamination_checks.R",
  "R/05_annotation_program_scoring.R",
  "R/06_DGE_cluster_markers_and_pseudobulk.R"
)

for (s in scripts) {
  message("Running: ", s)
  cmd <- sprintf("Rscript %s --config %s", shQuote(s), shQuote(args$config))
  status <- system(cmd)
  if (status != 0) stop("Pipeline failed at ", s, call. = FALSE)
}

message("Pipeline completed successfully.")
