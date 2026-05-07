#!/bin/bash -l
#SBATCH --job-name=retangio_02_qcfilter
#SBATCH --output=logs/02_qcfilter.out
#SBATCH --error=logs/02_qcfilter.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/02_qc_filter_cells.R
