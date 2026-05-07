#!/bin/bash -l
#SBATCH --job-name=retangio_06_sct
#SBATCH --output=logs/06_sct.out
#SBATCH --error=logs/06_sct.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/06_prepare_sct_objects.R
