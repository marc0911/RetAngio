#!/bin/bash -l
#SBATCH --job-name=retangio_09c_rpca
#SBATCH --output=logs/09c_finalize_ec_after_soupx_rpca.out
#SBATCH --error=logs/09c_finalize_ec_after_soupx_rpca.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs results/objects decisions

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/09c_finalize_ec_after_soupx_rpca.R
