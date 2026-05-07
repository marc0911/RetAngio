#!/bin/bash -l
#SBATCH --job-name=retangio_07_pca
#SBATCH --output=logs/07_pca.out
#SBATCH --error=logs/07_pca.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/07_run_pca_unintegrated.R
