# RetAngio reproducible Seurat v5 pipeline

This repository contains a configuration-driven, script-based pipeline for retinal endothelial-cell scRNA-seq re-analysis starting from an existing Seurat object.

## Scope

The pipeline performs:

1. Removal of non-EC contaminant clusters (`5,7,11,13` by default from `SCT_snn_res.0.4`).
2. Reclustering parameter sweep and selection.
3. Creation of a compact `BASE` Seurat snapshot.
4. QC and contamination/ambient RNA checks, with optional SoupX correction if raw droplet matrices are provided.
5. Program scoring-based cluster annotation using Zarkada and reference marker sets.
6. Marker discovery and pseudo-bulk differential expression.
7. PNG figure exports and structured outputs.

## Project structure

- `config/config.yml` — all runtime parameters.
- `R/00_utils.R` — logging, safe I/O, object checks, reproducibility helpers.
- `R/01_load_and_snapshot.R` — input loading, structure validation, object report.
- `R/02_remove_contaminants_and_recluster.R` — contaminant removal + parameter search.
- `R/03_build_BASE_object.R` — compact BASE object construction.
- `R/04_QC_and_contamination_checks.R` — QC + contamination checks + optional SoupX branch.
- `R/05_annotation_program_scoring.R` — marker/program scoring and annotation table.
- `R/06_DGE_cluster_markers_and_pseudobulk.R` — cluster markers + pseudo-bulk DGE.
- `R/run_all.R` — executes all modules in order.
- `Makefile` — convenience wrapper.

## Reproducibility and integrity controls

- Deterministic random seeds from config (`set.seed`).
- Parameterized outputs with explicit filenames.
- Session and package provenance export.
- Script-level logging and elapsed-time tracking.
- No hidden state: each step consumes prior explicit artifacts.

## Inputs required at runtime

Update `config/config.yml` with:

- `paths.input_seurat_rds` (required)
- `paths.zarkada_excel` (required for annotation script)
- Optional: `paths.raw_matrix_dirs` per sample for SoupX branch

## Running

```bash
make run
# or
Rscript R/run_all.R --config config/config.yml
```

Run individual steps if needed, e.g.:

```bash
Rscript R/02_remove_contaminants_and_recluster.R --config config/config.yml
```

## Package management (`renv`)

This repo includes bootstrap instructions in `scripts/setup_renv.R`.

```bash
Rscript scripts/setup_renv.R
```

If `renv.lock` is absent, this script initializes `renv` and snapshots package versions from the current environment. If present, it restores the lockfile.

## Key references

> Note: The Zarkada marker manuscript citation is intentionally left as a user-supplied config/docs item because multiple Zarkada vascular marker resources are used in practice and must be mapped to the exact Excel file delivered at runtime.

- Seurat v5 command and layer model documentation: Satija Lab Seurat reference docs (authoritative docs).
- Vanlandewijck et al. *Nature* (2018), PMID: 29443965, DOI: 10.1038/nature25739.
- Young & Behjati. SoupX. *GigaScience* (2020), DOI: 10.1093/gigascience/giaa151.
- Germain et al. scDblFinder. *F1000Research* (2021), DOI: 10.12688/f1000research.73600.2.
- Soneson & Robinson. Bias, robustness and scalability in single-cell DE methods. *Nat Methods* (2018), DOI: 10.1038/nmeth.4612.
- Crowell et al. muscat framework for multi-sample multi-group scRNA-seq analysis. *Nat Commun* (2020), DOI: 10.1038/s41467-020-19894-4.
- Chen et al. fgsea. *bioRxiv* preprint / package docs; pair with pathway databases used.

## Outputs

Outputs are written to `outputs/` (created automatically) and are not committed.
