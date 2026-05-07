#!/bin/bash -l
#SBATCH --job-name=retangio_08d_plot
#SBATCH --output=logs/08d_plot_ec_after_soupx_unintegrated_review.out
#SBATCH --error=logs/08d_plot_ec_after_soupx_unintegrated_review.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs plots/08d_ec_after_soupx_unintegrated_review

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/08d_plot_ec_after_soupx_unintegrated_review.R
