#!/bin/bash
#SBATCH --job-name=09s_ec_panel_v2
#SBATCH --output=logs/09s_plot_curated_overall_ec_panel_rpca_v2.%j.out
#SBATCH --error=logs/09s_plot_curated_overall_ec_panel_rpca_v2.%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/09s_plot_curated_overall_ec_panel_rpca_v2.R
