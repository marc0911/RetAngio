# RetAngio v4

RetAngio v4 is a reproducible scRNA-seq analysis repository for the mouse retinal OIR project. It is designed to start from raw Cell Ranger outputs, run cleanly on HPC with SLURM, and enforce explicit human decision checkpoints at key analysis stages.

## Current status

Only **Stage 0** is implemented at the moment:

- build a fresh merged Seurat object from sample-level filtered Cell Ranger matrices
- attach sample metadata
- save a resolved sample manifest
- export a per-sample build summary

No QC filtering, doublet removal, ambient RNA correction, clustering, integration, annotation, or DGE is implemented yet.

## Planned pipeline stages

0. Build fresh object from Cell Ranger outputs
1. Compute QC metrics
2. Filter low-quality cells
3. Detect/remove doublets
4. Ambient RNA correction
5. Unintegrated clustering
6. Remove non-endothelial cells
7. Compare integration methods: none vs RPCA vs Harmony
8. Choose final integration
9. Annotation / program scoring
10. DGE / pseudobulk

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
├── slurm/
├── data/
│   ├── raw_archives/
│   ├── extracted/
│   └── metadata/
├── results/
│   ├── objects/
│   ├── qc/
│   ├── doublets/
│   ├── ambient/
│   ├── unintegrated/
│   ├── integration_compare/
│   ├── annotation/
│   └── dge/
├── plots/
└── logs/
