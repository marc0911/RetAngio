#!/bin/bash
#SBATCH --job-name=09za_tip_ctrl_sd
#SBATCH --output=logs/09za_plot_tipcell_ctrl_only_sd_marker_panel_final.%j.out
#SBATCH --error=logs/09za_plot_tipcell_ctrl_only_sd_marker_panel_final.%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/09za_plot_tipcell_ctrl_only_sd_marker_panel_final.R
