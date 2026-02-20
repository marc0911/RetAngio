source("R/00_utils.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(edgeR)
  library(ggplot2)
})

# Method notes:
# - Cluster marker testing uses Seurat FindAllMarkers on RNA assay.
# - Pseudo-bulk DE is used to avoid pseudoreplication at cell level
#   (Soneson & Robinson 2018 DOI:10.1038/nmeth.4612; Crowell et al. 2020 DOI:10.1038/s41467-020-19894-4).
# - edgeR quasi-likelihood pipeline follows edgeR user guide best practices.

args <- parse_args()
cfg <- load_config(args$config)
set.seed(cfg$seed)
init_dirs(cfg)
script_name <- "06_DGE_cluster_markers_and_pseudobulk"

infile <- file.path(cfg$output$root, cfg$output$snapshots_dir, "05_annotated.rds")
if (!file.exists(infile)) stop("Missing 05 output.", call. = FALSE)
seu <- readRDS(infile)

cluster_col <- paste0("SCT_snn_res.", cfg$recluster$final_resolution)
if (!(cluster_col %in% colnames(seu@meta.data))) {
  cluster_col <- tail(grep("^SCT_snn_res\\.", colnames(seu@meta.data), value = TRUE), 1)
}
Idents(seu) <- seu@meta.data[[cluster_col]]

DefaultAssay(seu) <- "RNA"
marker_cfg <- cfg$dge$cluster_markers
markers <- FindAllMarkers(
  object = seu,
  assay = "RNA",
  slot = "data",
  only.pos = FALSE,
  test.use = marker_cfg$test_use,
  min.pct = marker_cfg$min_pct,
  logfc.threshold = marker_cfg$logfc_threshold,
  max.cells.per.ident = marker_cfg$max_cells_per_ident,
  verbose = FALSE
)
readr::write_csv(markers, file.path(cfg$output$root, cfg$output$tables_dir, "06_cluster_markers_up_down.csv"))

meta <- seu@meta.data %>% tibble::rownames_to_column("cell")
required_cols <- c(cluster_col, cfg$dge$pseudobulk$sample_column, cfg$dge$pseudobulk$condition_column)
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop(glue::glue("Missing required metadata for pseudo-bulk: {paste(missing_cols, collapse=', ')}"), call. = FALSE)
}

time_col <- cfg$dge$pseudobulk$timepoint_column
if (!(time_col %in% colnames(meta))) {
  meta[[time_col]] <- "NA"
}

counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
meta$sample_cluster <- paste(meta[[cfg$dge$pseudobulk$sample_column]], meta[[cluster_col]], sep = "__")

sample_cluster_info <- meta %>%
  count(sample_cluster, name = "n_cells") %>%
  left_join(meta %>% distinct(sample_cluster, .keep_all = TRUE), by = "sample_cluster") %>%
  filter(n_cells >= cfg$dge$pseudobulk$min_cells_per_sample_cluster)

keep_sc <- sample_cluster_info$sample_cluster
if (length(keep_sc) < 2) stop("Not enough pseudo-bulk replicates after min-cells filter.", call. = FALSE)

mat <- as.matrix(counts[, meta$sample_cluster %in% keep_sc])
sc_labels <- meta$sample_cluster[meta$sample_cluster %in% keep_sc]
agg <- rowsum(t(mat), group = sc_labels, reorder = FALSE)
agg <- t(agg)

pb_meta <- sample_cluster_info %>%
  select(sample_cluster,
         sample = all_of(cfg$dge$pseudobulk$sample_column),
         cluster = all_of(cluster_col),
         condition = all_of(cfg$dge$pseudobulk$condition_column),
         timepoint = all_of(time_col)) %>%
  distinct() %>%
  filter(sample_cluster %in% colnames(agg))

results_list <- list()
for (cl in unique(pb_meta$cluster)) {
  idx <- pb_meta$cluster == cl
  y <- DGEList(counts = agg[, idx, drop = FALSE], samples = pb_meta[idx, ])
  keep <- filterByExpr(y, group = y$samples$condition)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  design <- model.matrix(~ condition + timepoint, data = y$samples)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)

  coef_name <- grep("condition", colnames(design), value = TRUE)[1]
  if (!is.na(coef_name)) {
    qlf <- glmQLFTest(fit, coef = which(colnames(design) == coef_name))
    tab <- topTags(qlf, n = Inf)$table %>%
      tibble::rownames_to_column("gene") %>%
      mutate(cluster = cl, contrast = coef_name)
    results_list[[as.character(cl)]] <- tab

    volcano <- ggplot(tab, aes(x = logFC, y = -log10(PValue))) +
      geom_point(alpha = 0.4, size = 0.7) +
      geom_vline(xintercept = c(-1, 1), linetype = 2) +
      ggtitle(paste("Pseudo-bulk volcano cluster", cl))
    save_plot(volcano, paste0("06_pseudobulk_volcano_cluster_", cl, ".png"), cfg)
  }
}

pb_res <- bind_rows(results_list)
readr::write_csv(pb_res, file.path(cfg$output$root, cfg$output$tables_dir, "06_pseudobulk_edgeR_results.csv"))

sanity <- tibble::tibble(
  check = c("idents_match_final_cluster", "markers_nonempty", "pseudobulk_nonempty"),
  pass = c(
    identical(as.character(Idents(seu)), as.character(seu@meta.data[[cluster_col]])),
    nrow(markers) > 0,
    nrow(pb_res) > 0
  )
)
readr::write_csv(sanity, file.path(cfg$output$root, cfg$output$reports_dir, "06_sanity_checks.csv"))

write_session_info(cfg, script_name)
log_message("Completed DE analyses", cfg, script_name)
