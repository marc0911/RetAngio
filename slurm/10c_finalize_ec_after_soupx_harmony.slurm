#!/bin/bash -l
#SBATCH --job-name=retangio_10c_harmony
#SBATCH --output=logs/10c_finalize_ec_after_soupx_harmony.out
#SBATCH --error=logs/10c_finalize_ec_after_soupx_harmony.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs results/objects decisions

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/10c_finalize_ec_after_soupx_harmony.R
