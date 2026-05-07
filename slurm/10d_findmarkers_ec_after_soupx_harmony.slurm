#!/bin/bash -l
#SBATCH --job-name=retangio_10d_markers
#SBATCH --output=logs/10d_findmarkers_ec_after_soupx_harmony.out
#SBATCH --error=logs/10d_findmarkers_ec_after_soupx_harmony.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs results/markers

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/10d_findmarkers_ec_after_soupx_harmony.R
