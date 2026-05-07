#!/bin/bash
#SBATCH --job-name=09w_tip_overall_unscaled
#SBATCH --output=logs/09w_plot_tipcell_overall_curated_panel_unscaled.%j.out
#SBATCH --error=logs/09w_plot_tipcell_overall_curated_panel_unscaled.%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

# IMPORTANT:
# If your previous working overall tip-cell script used a different object path,
# replace ONLY the line below with that exact object path.
export INPUT_RDS="results/objects/09k_rpca_tipcell_resolution_sweep.rds"

Rscript scripts/09w_plot_tipcell_overall_curated_panel_unscaled.R
