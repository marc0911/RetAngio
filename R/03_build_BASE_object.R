source("R/00_utils.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# Method notes:
# - Keep RNA assay normalized data layer for DEG readiness per Seurat guidance.
# - JoinLayers is used for Seurat v5 layer coherence before DE-oriented workflows.

args <- parse_args()
cfg <- load_config(args$config)
set.seed(cfg$seed)
init_dirs(cfg)
script_name <- "03_build_BASE_object"

infile <- file.path(cfg$output$root, cfg$output$snapshots_dir, "02_post_filter_reclustered.rds")
if (!file.exists(infile)) stop("Missing 02 output.", call. = FALSE)
seu <- readRDS(infile)

meta_keep <- intersect(cfg$base_object$keep_meta_columns, colnames(seu@meta.data))
seu@meta.data <- seu@meta.data[, unique(c(meta_keep, colnames(seu@meta.data)[grepl("^SCT_snn_res\\.", colnames(seu@meta.data))])), drop = FALSE]

DefaultAssay(seu) <- "RNA"
seu <- JoinLayers(seu, assay = "RNA")
seu <- NormalizeData(seu, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)

seu <- compactify_seurat(seu, cfg)

base_path <- file.path(cfg$output$root, cfg$output$snapshots_dir, "03_BASE_object.rds")
checksum <- safe_save_rds(seu, base_path)
writeLines(checksum, file.path(cfg$output$root, cfg$output$reports_dir, "03_BASE_object.sha256.txt"))

has_data <- "data" %in% SeuratObject::Layers(seu[["RNA"]])
check_tbl <- tibble::tibble(
  check = c("rna_has_data_layer", "assays_subset"),
  pass = c(has_data, all(Assays(seu) %in% cfg$base_object$assays_to_keep))
)
readr::write_csv(check_tbl, file.path(cfg$output$root, cfg$output$reports_dir, "03_sanity_checks.csv"))

write_session_info(cfg, script_name)
log_message("Completed BASE object build", cfg, script_name)
