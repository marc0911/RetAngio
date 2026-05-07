#!/bin/bash -l
#SBATCH --job-name=retangio_07b_review
#SBATCH --output=logs/07b_prepare_cluster_review_after_soupx.out
#SBATCH --error=logs/07b_prepare_cluster_review_after_soupx.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs results/review

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/07b_prepare_cluster_review_after_soupx.R
