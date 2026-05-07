#!/bin/bash -l
#SBATCH --job-name=09z_ec_panel_ng
#SBATCH --output=logs/09z_plot_curated_overall_ec_panel_rpca_v2_nogrid.%j.out
#SBATCH --error=logs/09z_plot_curated_overall_ec_panel_rpca_v2_nogrid.%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/09z_plot_curated_overall_ec_panel_rpca_v2_nogrid.R
