#!/bin/bash -l
#SBATCH --job-name=retangio_03_doublets
#SBATCH --output=logs/03_doublets.out
#SBATCH --error=logs/03_doublets.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/03_detect_doublets.R
