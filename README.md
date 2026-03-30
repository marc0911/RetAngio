# RetAngio (v4 rebuild)

## Purpose
RetAngio is a reproducible, non-interactive scRNA-seq analysis repository for mouse retinal OIR data. This rebuild starts from raw Cell Ranger outputs and is designed for both local and HPC/SLURM execution.

## Project overview
This project currently tracks five samples:
- P7
- P12_CTRL
- P12_OIR
- P17_CTRL
- P17_OIR

The repository is organized as explicit stage scripts plus decision checkpoints.

## Current implementation status
Only **Stage 0** is currently implemented.

Implemented now:
- Build a fresh merged Seurat object from per-sample filtered matrices.
- Preserve sample-level provenance and metadata.
- Save required stage outputs to disk.

Not implemented yet:
- QC filtering
- Doublet detection/removal
- Ambient RNA correction
- Clustering/integration workflows
- Cell annotation/program scoring
- DGE/pseudobulk analysis

## Expected raw inputs
For each sample, extracted Cell Ranger archive content is expected to include:
- `raw_feature_bc_matrix` or `raw_feature_bc_matrix.h5`
- `sample_filtered_feature_bc_matrix` or `sample_filtered_feature_bc_matrix.h5`
- `raw_molecule_info.h5`

Recommended layout:
- Archives: `data/raw_archives/`
- Extracted sample folders/files: `data/extracted/`

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
```

## Configure `config/samples.csv`
1. Open `config/samples.csv`.
2. Keep one row per sample.
3. Set required stage-0 fields:
   - `sample_id`
   - `condition`
   - `timepoint`
   - `filtered_h5` (required for stage 0)
4. Fill `raw_h5` and `molecule_info_h5` if available (recommended for provenance).

For stage 0, only the filtered matrix input is required to run.

## Run locally
From the repository root:

```bash
Rscript R/00_build_fresh_object.R config/config.yml
```

If config path is omitted, default is `config/config.yml`.

## Run on cluster (SLURM)
Submit:

```bash
sbatch slurm/run_00_build_fresh_object.sh
```

Default SLURM resources are conservative (8 CPUs, 64G RAM, 12h) and may need adjustment.

## Stage-0 outputs
- `results/objects/00_raw_merged.rds`
- `results/qc/00_build_summary.csv`
- `data/metadata/samples_resolved.csv`

## Planned future stages
1. QC metric computation
2. QC filtering
3. Doublet detection/removal (scDblFinder)
4. Ambient RNA correction (SoupX)
5. Unintegrated clustering
6. Non-endothelial contamination removal
7. Integration comparison (none vs RPCA vs Harmony)
8. Integration selection
9. Annotation/program scoring
10. DGE/pseudobulk
