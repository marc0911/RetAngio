#!/bin/bash -l
#SBATCH --job-name=retangio_18_ec_markers
#SBATCH --output=logs/18_markers_ec_only_integrated_RNA.out
#SBATCH --error=logs/18_markers_ec_only_integrated_RNA.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo
mkdir -p logs results/markers

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/18_markers_ec_only_integrated_RNA.R
