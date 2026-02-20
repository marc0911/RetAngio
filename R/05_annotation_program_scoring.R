source("R/00_utils.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(purrr)
  library(stringr)
  library(tidyr)
})

# Method notes:
# - Marker-based vascular annotation references Vanlandewijck et al. Nature 2018.
#   PMID: 29443965; DOI: 10.1038/nature25739.
# - Module scoring defaults to AddModuleScore (Seurat docs).
# - If desired, users can switch to AUCell/UCell in future extensions.
# - IMPORTANT: Zarkada citation should be filled to match the exact marker Excel provided.

args <- parse_args()
cfg <- load_config(args$config)
set.seed(cfg$seed)
init_dirs(cfg)
script_name <- "05_annotation_program_scoring"

infile <- file.path(cfg$output$root, cfg$output$snapshots_dir, "04_qc_checked.rds")
if (!file.exists(infile)) stop("Missing 04 output.", call. = FALSE)
seu <- readRDS(infile)

zarkada_path <- cfg$paths$zarkada_excel
if (!file.exists(zarkada_path)) stop("Zarkada Excel file missing.", call. = FALSE)

sheet_names <- readxl::excel_sheets(zarkada_path)
zarkada_tbl <- purrr::map_dfr(sheet_names, function(sh) {
  readxl::read_excel(zarkada_path, sheet = sh) %>%
    janitor::clean_names() %>%
    mutate(source_sheet = sh)
}, .id = "sheet_id")

candidate_gene_col <- names(zarkada_tbl)[str_detect(names(zarkada_tbl), "gene|symbol")][1]
candidate_set_col <- names(zarkada_tbl)[str_detect(names(zarkada_tbl), "set|program|type|class|cluster")][1]
if (is.na(candidate_gene_col) || is.na(candidate_set_col)) {
  stop("Could not infer gene/program columns from Zarkada Excel. Please standardize columns.", call. = FALSE)
}

zarkada_sets <- zarkada_tbl %>%
  transmute(program = as.character(.data[[candidate_set_col]]), gene = as.character(.data[[candidate_gene_col]])) %>%
  filter(!is.na(program), !is.na(gene), gene != "") %>%
  distinct() %>%
  group_by(program) %>%
  summarise(genes = list(unique(gene)), n_genes = n(), .groups = "drop") %>%
  filter(n_genes >= cfg$annotation$min_geneset_size)

vanlandewijck_sets <- tibble::tribble(
  ~program, ~gene,
  "arterial_reference", "Gja5",
  "arterial_reference", "Fbln5",
  "arterial_reference", "Efnb2",
  "venous_reference", "Nr2f2",
  "venous_reference", "Vwf",
  "venous_reference", "Selp",
  "capillary_reference", "Kdr",
  "capillary_reference", "Cd36",
  "capillary_reference", "Rgcc"
) %>% group_by(program) %>% summarise(genes = list(unique(gene)), n_genes = n(), .groups = "drop")

all_sets <- bind_rows(zarkada_sets, vanlandewijck_sets) %>%
  mutate(genes = purrr::map(genes, ~ intersect(.x, rownames(seu)))) %>%
  mutate(n_present = purrr::map_int(genes, length)) %>%
  filter(n_present >= cfg$annotation$min_geneset_size)

if (nrow(all_sets) == 0) stop("No gene sets available after filtering/presence checks.", call. = FALSE)

seu <- AddModuleScore(
  object = seu,
  features = all_sets$genes,
  name = "ProgScore",
  assay = "RNA",
  search = TRUE
)

score_cols <- grep("^ProgScore", colnames(seu@meta.data), value = TRUE)
score_map <- tibble::tibble(program = all_sets$program, score_col = score_cols[seq_len(nrow(all_sets))])

cluster_col <- paste0("SCT_snn_res.", cfg$recluster$final_resolution)
if (!(cluster_col %in% colnames(seu@meta.data))) {
  cluster_col <- tail(grep("^SCT_snn_res\\.", colnames(seu@meta.data), value = TRUE), 1)
}

cluster_scores <- seu@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  select(all_of(c(cluster_col, score_map$score_col))) %>%
  pivot_longer(cols = all_of(score_map$score_col), names_to = "score_col", values_to = "score") %>%
  left_join(score_map, by = "score_col") %>%
  group_by(.data[[cluster_col]], program) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

annotations <- cluster_scores %>%
  group_by(.data[[cluster_col]]) %>%
  slice_max(order_by = mean_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(cluster = all_of(cluster_col), top_program = program, top_program_score = mean_score)

readr::write_csv(cluster_scores, file.path(cfg$output$root, cfg$output$tables_dir, "05_cluster_program_scores.csv"))
readr::write_csv(annotations, file.path(cfg$output$root, cfg$output$tables_dir, "05_cluster_annotations.csv"))

for (sc in head(score_map$score_col, 12)) {
  p <- FeaturePlot(seu, features = sc, reduction = "umap") + ggplot2::ggtitle(sc)
  save_plot(p, paste0("05_featureplot_", sc, ".png"), cfg)
}

out <- file.path(cfg$output$root, cfg$output$snapshots_dir, "05_annotated.rds")
safe_save_rds(seu, out)
write_session_info(cfg, script_name)
log_message("Completed annotation + program scoring", cfg, script_name)
