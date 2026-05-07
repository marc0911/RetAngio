#!/bin/bash -l
#SBATCH --job-name=retangio_08f_markers
#SBATCH --output=logs/08f_findmarkers_ec_after_soupx_unintegrated.out
#SBATCH --error=logs/08f_findmarkers_ec_after_soupx_unintegrated.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs results/markers

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/08f_findmarkers_ec_after_soupx_unintegrated.R
