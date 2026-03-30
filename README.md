# RetAngio (v4 rebuild)

## Purpose
RetAngio is a reproducible, non-interactive scRNA-seq analysis repository for mouse retinal OIR data. This rebuild starts from raw Cell Ranger outputs and is designed for local execution and HPC/SLURM environments.

## Project overview
This project currently tracks five samples:
- P7
- P12_CTRL
- P12_OIR
- P17_CTRL
- P17_OIR

The repository is designed around explicit stage boundaries and YAML decision checkpoints.

## Expected raw inputs
For each sample, extracted Cell Ranger archive content should include:
- `raw_feature_bc_matrix` or `raw_feature_bc_matrix.h5`
- `sample_filtered_feature_bc_matrix` or `sample_filtered_feature_bc_matrix.h5`
- `raw_molecule_info.h5`

Place raw archives under `data/raw_archives/` and extracted per-sample files under `data/extracted/`.

## Implemented stage (current)
Only **Stage 0** is implemented now:
- Build a fresh merged Seurat object from per-sample filtered matrices.
- Preserve sample metadata/provenance.
- Save merged object and summary tables to disk.

No QC filtering, doublet detection, ambient RNA correction, clustering, integration, annotation, or DGE is performed in stage 0.

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

## Configure sample manifest (`config/samples.csv`)
1. Open `config/samples.csv`.
2. Confirm `sample_id`, `condition`, and `timepoint` values.
3. Replace placeholder file paths for:
   - `filtered_h5`
   - `raw_h5`
   - `molecule_info_h5`
4. Ensure every path exists before running stage 0.

The manifest is treated as the source of truth.

## Run locally
From repository root:

```bash
Rscript R/00_build_fresh_object.R config/config.yml
```

If config path is omitted, the script defaults to `config/config.yml`.

## Run on cluster (SLURM)
Submit the stage-0 wrapper:

```bash
sbatch slurm/run_00_build_fresh_object.sh
```

The script requests conservative defaults (8 CPUs, 64G RAM, 12h). Tune these for your system.

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
9. Annotation and program scoring
10. DGE and pseudobulk analysis
