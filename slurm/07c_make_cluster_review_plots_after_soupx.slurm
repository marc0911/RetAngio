#!/bin/bash -l
#SBATCH --job-name=retangio_07c_plots
#SBATCH --output=logs/07c_make_cluster_review_plots_after_soupx.out
#SBATCH --error=logs/07c_make_cluster_review_plots_after_soupx.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs plots/07c_cluster_review_after_soupx

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/07c_make_cluster_review_plots_after_soupx.R
