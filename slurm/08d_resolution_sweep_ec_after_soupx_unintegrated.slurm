#!/bin/bash -l
#SBATCH --job-name=retangio_08d_res
#SBATCH --output=logs/08d_resolution_sweep_ec_after_soupx_unintegrated.out
#SBATCH --error=logs/08d_resolution_sweep_ec_after_soupx_unintegrated.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs plots/08d_ec_after_soupx_unintegrated_resolution_sweep results/objects

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/08d_resolution_sweep_ec_after_soupx_unintegrated.R
