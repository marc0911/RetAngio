#!/bin/bash -l
#SBATCH --job-name=retangio_10a_harmony
#SBATCH --output=logs/10a_integrate_ec_after_soupx_harmony.out
#SBATCH --error=logs/10a_integrate_ec_after_soupx_harmony.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs results/objects

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/10a_integrate_ec_after_soupx_harmony.R
