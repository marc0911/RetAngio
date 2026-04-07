#!/bin/bash -l
#SBATCH --job-name=retangio_11_subset_ec
#SBATCH --output=logs/11_subset_ec.out
#SBATCH --error=logs/11_subset_ec.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo
mkdir -p logs results/objects results/metadata

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/11_subset_ec_from_unintegrated_clusters.R
