#!/bin/bash -l
#SBATCH --job-name=retangio_08d_legends
#SBATCH --output=logs/08d_regenerate_legends_only.out
#SBATCH --error=logs/08d_regenerate_legends_only.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/08d_regenerate_legends_only.R
