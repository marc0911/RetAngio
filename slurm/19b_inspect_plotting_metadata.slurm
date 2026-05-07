#!/bin/bash -l
#SBATCH --job-name=retangio_plot_meta
#SBATCH --output=logs/19b_inspect_plotting_metadata.out
#SBATCH --error=logs/19b_inspect_plotting_metadata.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/19b_inspect_plotting_metadata.R
