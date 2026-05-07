#!/bin/bash -l
#SBATCH --job-name=retangio_05b_soupx
#SBATCH --output=logs/05b_run_soupx.out
#SBATCH --error=logs/05b_run_soupx.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs results/objects results/ambient

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/05b_run_soupx.R
