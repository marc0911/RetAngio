#!/bin/bash -l
#SBATCH --job-name=retangio_08a_sct
#SBATCH --output=logs/08a_prepare_sct_ec_after_soupx.out
#SBATCH --error=logs/08a_prepare_sct_ec_after_soupx.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs results/objects

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/08a_prepare_sct_ec_after_soupx.R
