#!/bin/bash -l
#SBATCH --job-name=retangio_install_readxl
#SBATCH --output=logs/00_install_readxl.out
#SBATCH --error=logs/00_install_readxl.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/00_install_readxl.R
