source("R/00_utils.R")
suppressPackageStartupMessages({
  library(readr)
})

args <- parse_args()
cfg <- load_config(args$config)
set.seed(cfg$seed)
init_dirs(cfg)
script_name <- "01_load_and_snapshot"

log_message("Starting load + snapshot step", cfg, script_name)

input_rds <- cfg$paths$input_seurat_rds
if (!file.exists(input_rds)) {
  stop(glue::glue("Input Seurat RDS not found: {input_rds}"), call. = FALSE)
}

seu <- time_it(readRDS(input_rds), "Read input Seurat object", cfg, script_name)
validate_seurat(seu, cfg)

report <- object_size_report(seu)
report_path <- file.path(cfg$output$root, cfg$output$reports_dir, "01_input_object_size_report.csv")
readr::write_csv(report, report_path)

snapshot_path <- file.path(cfg$output$root, cfg$output$snapshots_dir, "01_input_validated.rds")
checksum <- safe_save_rds(seu, snapshot_path)

checksum_path <- file.path(cfg$output$root, cfg$output$reports_dir, "01_input_validated.sha256.txt")
writeLines(checksum, checksum_path)

write_session_info(cfg, script_name)
log_message("Completed load + snapshot step", cfg, script_name)
