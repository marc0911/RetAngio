#!/bin/bash -l
#SBATCH --job-name=retangio_14_ec_umap
#SBATCH --output=logs/14_cluster_umap_ec_only_unintegrated.out
#SBATCH --error=logs/14_cluster_umap_ec_only_unintegrated.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo
mkdir -p logs results/objects results/clusters plots/umap

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/14_cluster_umap_ec_only_unintegrated.R
