#!/bin/bash -l
#SBATCH --job-name=retangio_07a_markers
#SBATCH --output=logs/07a_findmarkers_unintegrated_after_soupx.out
#SBATCH --error=logs/07a_findmarkers_unintegrated_after_soupx.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs results/markers

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/07a_findmarkers_unintegrated_after_soupx.R
