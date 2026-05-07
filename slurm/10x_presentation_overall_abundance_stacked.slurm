#!/bin/bash -l
#SBATCH --job-name=10x_overall_abund
#SBATCH --output=logs/10x_presentation_overall_abundance_stacked.%j.out
#SBATCH --error=logs/10x_presentation_overall_abundance_stacked.%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs:${R_LIBS_USER:-}

Rscript scripts/10x_presentation_overall_abundance_stacked.R
