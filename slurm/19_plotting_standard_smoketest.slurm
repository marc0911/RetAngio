#!/bin/bash -l
#SBATCH --job-name=retangio_plot_std
#SBATCH --output=logs/19_plotting_standard_smoketest.out
#SBATCH --error=logs/19_plotting_standard_smoketest.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs plots/19_plotting_standard_smoketest

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/19_plotting_standard_smoketest.R
