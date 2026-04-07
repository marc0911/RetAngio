# RetAngio v4

RetAngio v4 is a reproducible scRNA-seq analysis repository for the mouse retinal OIR project. It is designed to start from raw Cell Ranger outputs, run cleanly on HPC with SLURM, and enforce explicit human decision checkpoints at key analysis stages.

## Current status

The repository now contains a working Seurat v5 / SLURM pipeline implemented through the current **EC-only branch of analysis**.

### Implemented stages so far

- **Stage 0**: build fresh merged Seurat object from Cell Ranger filtered matrices
- **Stage 1**: compute QC metrics
- **Stage 2**: filter low-quality cells
- **Stage 3**: detect doublets
- **Stage 4**: remove doublets
- **Stage 5**: post-doublet QC check
- **Stage 6**: prepare per-sample SCT objects
- **Stage 7**: run unintegrated PCA
- **Stage 8**: unintegrated clustering and UMAP
- **Stage 9**: RPCA integration fast-track
- **Stage 10**: contaminant review and visualization
- **Stage 11**: subset EC-enriched cells from the unintegrated object
- **Stage 12**: prepare EC-only SCT objects
- **Stage 13**: run EC-only unintegrated PCA
- **Stage 14**: EC-only unintegrated clustering and UMAP
- **Stage 15**: EC-only marker calculation on RNA assay
- **Stage 16**: EC-only RPCA integration
- **Stage 17**: EC-only integrated clustering and UMAP
- **Stage 17b**: replot integrated UMAPs with harmonized styling
- **Stage 18**: EC-only integrated marker calculation on RNA assay

## Current analysis logic

The present repository reflects a two-step strategy:

1. perform initial clustering on the full filtered dataset,
2. identify and exclude non-endothelial contaminants,
3. rebuild an **EC-enriched** object,
4. rerun normalization, PCA, clustering, and integration on the EC-only branch.

This is the currently active and preferred analysis branch.

## Current key decisions

The following decision files are currently part of the workflow:

- `decisions/01_qc_thresholds.yml`
- `decisions/05_unintegrated_dims_resolution.yml`
- `decisions/06_contaminant_clusters.yml`
- `decisions/08_integration_choice.yml`

These files document the major human checkpoints used to guide the pipeline.

## Expected raw inputs

For each sample, the project expects Cell Ranger-derived files such as:

- `sample_filtered_feature_bc_matrix.h5`
- `raw_feature_bc_matrix.h5`
- `raw_molecule_info.h5`

The sample manifest lives in `config/samples.csv`.

## Repository structure

```text
RetAngio/
├── config/
├── decisions/
├── R/
├── scripts/
├── slurm/
├── data/
│   ├── raw_archives/
│   ├── extracted/
│   └── metadata/
├── results/
├── plots/
├── logs/
└── README.md
