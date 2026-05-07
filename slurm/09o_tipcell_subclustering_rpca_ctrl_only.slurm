#!/bin/bash -l
#SBATCH --job-name=09o_tip_ctrl
#SBATCH --output=logs/09o_tipcell_subclustering_rpca_ctrl_only_%j.out
#SBATCH --error=logs/09o_tipcell_subclustering_rpca_ctrl_only_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/09o_tipcell_subclustering_rpca_ctrl_only.R
