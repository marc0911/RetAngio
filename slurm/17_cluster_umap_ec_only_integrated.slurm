#!/bin/bash -l
#SBATCH --job-name=retangio_17_ec_int_umap
#SBATCH --output=logs/17_cluster_umap_ec_only_integrated.out
#SBATCH --error=logs/17_cluster_umap_ec_only_integrated.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo
mkdir -p logs results/objects results/clusters plots/integration

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/17_cluster_umap_ec_only_integrated.R
