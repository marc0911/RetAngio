#!/bin/bash -l
#SBATCH --job-name=retangio_10i_plot
#SBATCH --output=logs/10i_plot_annotated_ec_after_soupx_harmony.out
#SBATCH --error=logs/10i_plot_annotated_ec_after_soupx_harmony.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs plots/10i_ec_after_soupx_harmony_annotated

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/10i_plot_annotated_ec_after_soupx_harmony.R
