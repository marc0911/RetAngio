#!/bin/bash -l
#SBATCH --job-name=retangio_00_build
#SBATCH --output=logs/00_build.out
#SBATCH --error=logs/00_build.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/00_build_fresh_object.R config/config.yml
