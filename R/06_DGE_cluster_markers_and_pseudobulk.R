source("R/00_utils.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(edgeR)
  library(ggplot2)
  library(Matrix)
  library(tibble)
})

aggregate_cluster_sample <- function(counts, meta, cluster_id, sample_col) {
  idx <- meta$final_cluster == cluster_id
  sub_counts <- counts[, idx, drop = FALSE]
  sub_meta <- meta[idx, , drop = FALSE]
  mm <- Matrix::sparse.model.matrix(~ 0 + sub_meta[[sample_col]])
  colnames(mm) <- gsub("sub_meta\\[\\[sample_col\\]\\]", "", colnames(mm))
  agg <- sub_counts %*% mm
  list(counts = agg, samples = colnames(agg))
}

plot_volcano <- function(tbl, out_png) {
  p <- ggplot(tbl, aes(logFC, -log10(PValue))) +
    geom_point(alpha = 0.4, size = 0.8) +
    theme_bw()
  ggsave(out_png, p, width = 6, height = 5, dpi = 300)
}

plot_ma <- function(tbl, out_png) {
  p <- ggplot(tbl, aes(logCPM, logFC)) +
    geom_point(alpha = 0.4, size = 0.8) +
    theme_bw()
  ggsave(out_png, p, width = 6, height = 5, dpi = 300)
}

main <- function() {
  cfg <- read_config()
  outdir <- cfg$project$output_dir
  ensure_dir(outdir)

  seu <- load_object(file.path(outdir, "05_annotated"), cfg)
  check_required_layers(seu, "RNA")

  sample_col <- cfg$input$sample_column
  condition_col <- cfg$input$condition_column
  timepoint_col <- cfg$input$timepoint_column

  needed <- c(sample_col, condition_col, "final_cluster")
  if (!all(needed %in% colnames(seu@meta.data))) {
    fail_fast("Missing required metadata columns: {paste(setdiff(needed, colnames(seu@meta.data)), collapse=', ')}")
  }

  # Cluster markers on RNA assay for interpretability (Seurat DE docs)
  DefaultAssay(seu) <- "RNA"
  Idents(seu) <- "final_cluster"
  markers <- FindAllMarkers(
    seu,
    assay = "RNA",
    only.pos = FALSE,
    test.use = cfg$markers$test_use,
    logfc.threshold = cfg$markers$logfc_threshold,
    min.pct = cfg$markers$min_pct,
    max.cells.per.ident = cfg$markers$max_cells_per_ident,
    verbose = FALSE
  ) |>
    mutate(p_val_adj = p.adjust(p_val, method = "BH"))
  write_csv(markers, file.path(outdir, "06_cluster_markers_all.csv"))
  write_csv(filter(markers, avg_log2FC > 0), file.path(outdir, "06_cluster_markers_up.csv"))
  write_csv(filter(markers, avg_log2FC < 0), file.path(outdir, "06_cluster_markers_down.csv"))

  counts <- get_layer_data(seu, "RNA", "counts", fallback_slot = "counts")
  meta <- seu@meta.data |>
    mutate(cell = rownames(seu@meta.data)) |>
    as_tibble()

  sample_condition <- meta |>
    distinct(.data[[sample_col]], .data[[condition_col]], .keep_all = TRUE)
  if (nrow(sample_condition) < dplyr::n_distinct(meta[[sample_col]])) {
    fail_fast("sample_column maps to multiple condition labels; sample_column must be biological replicate ID")
  }

  rep_counts <- sample_condition |>
    count(.data[[condition_col]], name = "n_reps")
  write_csv(rep_counts, file.path(outdir, "06_replicates_per_condition.csv"))

  if (any(rep_counts$n_reps < cfg$pseudobulk$min_replicates_per_condition)) {
    warn <- "Insufficient biological replicates for inferential pseudobulk testing; exporting descriptive fold-change only"
    warning(warn, call. = FALSE)
  }

  clusters <- sort(unique(meta$final_cluster))
  qc_rows <- list()

  for (cl in clusters) {
    agg <- aggregate_cluster_sample(counts, meta, cl, sample_col)
    agg_counts <- agg$counts
    smeta <- sample_condition |>
      filter(.data[[sample_col]] %in% agg$samples)

    smeta <- smeta[match(colnames(agg_counts), smeta[[sample_col]]), , drop = FALSE]
    if (!identical(colnames(agg_counts), smeta[[sample_col]])) {
      fail_fast("Sample alignment failed for cluster {cl}")
    }

    smeta[[condition_col]] <- factor(smeta[[condition_col]])
    if (!cfg$input$baseline_condition %in% levels(smeta[[condition_col]])) {
      qc_rows[[length(qc_rows) + 1]] <- tibble(cluster = cl, status = "skipped", reason = "baseline condition absent")
      next
    }
    smeta[[condition_col]] <- relevel(smeta[[condition_col]], ref = cfg$input$baseline_condition)

    use_timepoint <- FALSE
    if (timepoint_col %in% colnames(smeta)) {
      tp <- factor(smeta[[timepoint_col]])
      confounded <- length(unique(interaction(tp, smeta[[condition_col]]))) == nlevels(tp) + nlevels(smeta[[condition_col]]) - 1
      use_timepoint <- nlevels(tp) > 1 && !confounded
      smeta[[timepoint_col]] <- tp
    }

    design <- if (use_timepoint) {
      model.matrix(as.formula(paste("~", condition_col, "+", timepoint_col)), data = smeta)
    } else {
      model.matrix(as.formula(paste("~", condition_col)), data = smeta)
    }

    if (qr(design)$rank < ncol(design)) {
      qc_rows[[length(qc_rows) + 1]] <- tibble(cluster = cl, status = "skipped", reason = "rank-deficient design")
      next
    }

    if (any(table(smeta[[condition_col]]) < cfg$pseudobulk$min_replicates_per_condition)) {
      avg_by_cond <- sapply(levels(smeta[[condition_col]]), function(cd) {
        Matrix::rowMeans(agg_counts[, smeta[[condition_col]] == cd, drop = FALSE])
      })
      if (is.null(dim(avg_by_cond))) {
        qc_rows[[length(qc_rows) + 1]] <- tibble(cluster = cl, status = "skipped", reason = "insufficient replicates")
        next
      }
      fc <- log2((avg_by_cond[, 2] + 1) / (avg_by_cond[, 1] + 1))
      write_csv(tibble(gene = rownames(agg_counts), log2FC_descriptive = fc), file.path(outdir, paste0("06_pseudobulk_cluster_", cl, "_descriptive_fc.csv")))
      qc_rows[[length(qc_rows) + 1]] <- tibble(cluster = cl, status = "skipped", reason = "insufficient replicates")
      next
    }

    y <- DGEList(counts = agg_counts, samples = smeta)
    keep <- filterByExpr(y, design = design)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design = design)
    fit <- glmQLFit(y, design, robust = TRUE)

    coef_name <- grep(paste0("^", condition_col), colnames(design), value = TRUE)
    coef_name <- coef_name[coef_name != condition_col]
    if (length(coef_name) != 1) {
      qc_rows[[length(qc_rows) + 1]] <- tibble(cluster = cl, status = "skipped", reason = "condition coefficient ambiguous")
      next
    }

    qlf <- glmQLFTest(fit, coef = which(colnames(design) == coef_name))
    tab <- topTags(qlf, n = Inf)$table |>
      rownames_to_column("gene")

    write_csv(tab, file.path(outdir, paste0("06_pseudobulk_cluster_", cl, "_full.csv")))
    write_csv(head(tab, cfg$pseudobulk$top_n), file.path(outdir, paste0("06_pseudobulk_cluster_", cl, "_top", cfg$pseudobulk$top_n, ".csv")))
    plot_volcano(tab, file.path(outdir, paste0("06_pseudobulk_cluster_", cl, "_volcano.png")))
    plot_ma(tab, file.path(outdir, paste0("06_pseudobulk_cluster_", cl, "_MA.png")))

    qc_rows[[length(qc_rows) + 1]] <- tibble(cluster = cl, status = "tested", reason = NA_character_)
  }

  qc_tbl <- bind_rows(qc_rows)
  write_csv(qc_tbl, file.path(outdir, "06_pseudobulk_qc_table.csv"))

  if (!any(qc_tbl$status == "tested") && !all(qc_tbl$status == "skipped")) {
    fail_fast("No tested clusters and no explicit skipped reasons recorded")
  }

  # Method references:
  # - Pseudobulk rationale: Soneson & Robinson 2018 (PMID: 29481549)
  # - Multi-sample multi-group DS context: Crowell et al. 2020 muscat (PMID: 32516394)
  # - edgeR QL workflow: edgeR User's Guide (Bioconductor authoritative documentation)
  log_info("06_DGE_cluster_markers_and_pseudobulk complete")
}

if (sys.nframe() == 0) {
  main()
}
