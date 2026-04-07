#!/bin/bash -l
#SBATCH --job-name=retangio_15_ec_markers
#SBATCH --output=logs/15_findmarkers_ec_only_rna.out
#SBATCH --error=logs/15_findmarkers_ec_only_rna.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo
mkdir -p logs results/objects results/markers

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/15_findmarkers_ec_only_rna.R
