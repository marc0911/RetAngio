#!/bin/bash -l
#SBATCH --job-name=09k_tip_rpca
#SBATCH --output=logs/09k_tipcell_subclustering_rpca_%j.out
#SBATCH --error=logs/09k_tipcell_subclustering_rpca_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/09k_tipcell_subclustering_rpca.R
