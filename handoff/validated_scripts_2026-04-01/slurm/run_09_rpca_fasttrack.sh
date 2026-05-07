#!/bin/bash -l
#SBATCH --job-name=retangio_09_rpca
#SBATCH --output=logs/09_rpca.out
#SBATCH --error=logs/09_rpca.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/09_rpca_fasttrack.R
