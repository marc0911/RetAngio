#!/bin/bash -l
#SBATCH --job-name=retangio_08b_pca
#SBATCH --output=logs/08b_run_pca_ec_after_soupx_unintegrated.out
#SBATCH --error=logs/08b_run_pca_ec_after_soupx_unintegrated.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs results/objects

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/08b_run_pca_ec_after_soupx_unintegrated.R
