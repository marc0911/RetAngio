#!/bin/bash -l
#SBATCH --job-name=tippres_overall
#SBATCH --output=logs/10x_presentation_plots_tipcells_overall.%j.out
#SBATCH --error=logs/10x_presentation_plots_tipcells_overall.%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/10x_presentation_plots_tipcells_overall.R
