#!/bin/bash -l
#SBATCH --job-name=09zc_ctrl_tip_fix
#SBATCH --output=logs/09zc_fix_ctrl_tip_sd_dotplot_labels_only.%j.out
#SBATCH --error=logs/09zc_fix_ctrl_tip_sd_dotplot_labels_only.%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/09zc_fix_ctrl_tip_sd_dotplot_labels_only.R
