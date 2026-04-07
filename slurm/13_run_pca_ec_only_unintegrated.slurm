#!/bin/bash -l
#SBATCH --job-name=retangio_13_ec_pca
#SBATCH --output=logs/13_run_pca_ec_only_unintegrated.out
#SBATCH --error=logs/13_run_pca_ec_only_unintegrated.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo
mkdir -p logs results/objects results/pca plots/pca

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/13_run_pca_ec_only_unintegrated.R
