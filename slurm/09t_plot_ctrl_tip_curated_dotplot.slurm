#!/bin/bash -l
#SBATCH --job-name=09t_ctrl_tip_dot
#SBATCH --output=logs/09t_plot_ctrl_tip_curated_dotplot.%j.out
#SBATCH --error=logs/09t_plot_ctrl_tip_curated_dotplot.%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/09t_plot_ctrl_tip_curated_dotplot.R
