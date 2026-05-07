#!/bin/bash -l
#SBATCH --job-name=retangio_08e_final
#SBATCH --output=logs/08e_finalize_ec_after_soupx_unintegrated.out
#SBATCH --error=logs/08e_finalize_ec_after_soupx_unintegrated.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs results/objects

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/08e_finalize_ec_after_soupx_unintegrated.R
