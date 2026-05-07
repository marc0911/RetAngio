#!/bin/bash -l
#SBATCH --job-name=09zb_tip_overall_curated
#SBATCH --output=logs/09zb_plot_tipcell_overall_sd_and_curated_pathological_panels.%j.out
#SBATCH --error=logs/09zb_plot_tipcell_overall_sd_and_curated_pathological_panels.%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/09zb_plot_tipcell_overall_sd_and_curated_pathological_panels.R
