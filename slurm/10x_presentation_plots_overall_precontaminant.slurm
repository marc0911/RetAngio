#!/bin/bash -l
#SBATCH --job-name=10x_pres_overall
#SBATCH --output=logs/10x_presentation_plots_overall_precontaminant.%j.out
#SBATCH --error=logs/10x_presentation_plots_overall_precontaminant.%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs:${R_LIBS_USER:-}

Rscript scripts/10x_presentation_plots_overall_precontaminant.R
