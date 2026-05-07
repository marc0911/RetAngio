#!/bin/bash -l
#SBATCH --job-name=retangio_10b_harmony
#SBATCH --output=logs/10b_resolution_sweep_ec_after_soupx_harmony.out
#SBATCH --error=logs/10b_resolution_sweep_ec_after_soupx_harmony.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs results/objects plots/10b_ec_after_soupx_harmony_resolution_sweep

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/10b_resolution_sweep_ec_after_soupx_harmony.R
