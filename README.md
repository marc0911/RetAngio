# RetAngio v2 pipeline

Reproducible Seurat v5.3.0 workflow for retinal endothelial-cell scRNA-seq with strict pseudobulk inference checks.

## Repository structure

- `config/config.yml`: single source of configuration.
- `R/00_utils.R`: shared helpers (logging, fail-fast checks, layer checks, IO abstraction, checksums).
- `R/01_load_and_snapshot.R`: input object loading and reproducibility snapshots.
- `R/02_remove_contaminants_and_reintegrate.R`: contaminant removal + SCT/RPCA reintegration.
- `R/03_build_BASE_object.R`: minimal BASE object snapshot.
- `R/04_QC_and_contamination_checks.R`: QC and contamination score summaries.
- `R/05_annotation_program_scoring.R`: program scoring and annotation summary.
- `R/06_DGE_cluster_markers_and_pseudobulk.R`: markers + replicate-aware pseudobulk DGE.
- `R/run_all.R`: end-to-end runner with terminal safety checks.
- `outputs/`: generated artifacts.

## Quick start

```bash
make run
```

## Serialization backends (portable by default)

Configured in `config/config.yml` under `io`:

- `rds` (default): `saveRDS()` / `readRDS()` from base R (portable, no optional packages).
- `qs2` (optional): `qs2::qs_save()` / `qs2::qs_read()` when `qs2` is installed.
- `fst_df` (optional): `fst::write_fst()` / `fst::read_fst()` for **data frames only**.

Behavioral guarantees:

- Seurat objects always save safely; if `fst_df` is configured for a non-data.frame object, pipeline falls back to RDS with warning.
- If `qs2`/`fst` is configured but package is unavailable, pipeline falls back to RDS with warning.
- Output extensions are backend-consistent: `.rds`, `.qs2`, `.fst`.

## Statistical validity guardrails

This pipeline enforces biological replication requirements before pseudobulk inference:

1. `sample_column` must represent biological replicate (not just condition label).
2. Minimum replicates per condition (`pseudobulk.min_replicates_per_condition`, default `2`) are required.
3. If replication is insufficient, inferential edgeR testing is skipped and descriptive log2 fold-change files are exported only.
4. Per-cluster tested/skipped status and reasons are written to `outputs/06_pseudobulk_qc_table.csv`.

## Reproducibility outputs

- `outputs/config_used.yml`
- `outputs/sessionInfo.txt`
- `outputs/package_versions.csv`
- `outputs/03_BASE_object_checksum.csv`
- `outputs/00_io_sanity_checks.csv`
- `renv.lock`

## Method references

- Seurat v5 integration docs (SCT + layer-based integration, RPCA):
  - https://satijalab.org/seurat/articles/seurat5_integration.html
- Seurat differential expression docs:
  - https://satijalab.org/seurat/reference/findallmarkers
- Pseudobulk rationale:
  - Soneson & Robinson 2018. PMID: 29481549
- Multi-sample DS analysis context:
  - Crowell et al. 2020 (muscat). PMID: 32516394
- SoupX ambient RNA correction:
  - Young & Behjati 2020. PMID: 33367645
- UCell rank-based scoring:
  - Andreatta & Carmona 2021. PMID: 34285779
- edgeR quasi-likelihood workflow (authoritative documentation):
  - edgeR User's Guide (Bioconductor): https://bioconductor.org/packages/edgeR
- Base R serialization docs:
  - saveRDS/readRDS: https://stat.ethz.ch/R-manual/R-devel/library/base/html/readRDS.html
- qs2 documentation (CRAN index + vignettes/manual):
  - https://cran.r-project.org/package=qs2
- fst documentation (official package site):
  - https://www.fstpackage.org/

## Notes

- Do not rerun `SCTransform` after contaminant removal unless raw counts changed.
- For Seurat v5, explicit layer usage is enforced (`counts`, `data`); legacy slot fallback is warning-gated.
- Scripts are written to run independently from clean sessions.
