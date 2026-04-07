#!/bin/bash -l
#SBATCH --job-name=retangio_17b_replot
#SBATCH --output=logs/17b_replot_ec_only_integrated.out
#SBATCH --error=logs/17b_replot_ec_only_integrated.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo
mkdir -p logs plots/integration

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/17b_replot_ec_only_integrated.R
