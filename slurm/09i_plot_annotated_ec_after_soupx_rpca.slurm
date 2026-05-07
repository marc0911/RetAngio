#!/bin/bash -l
#SBATCH --job-name=retangio_09i_plot
#SBATCH --output=logs/09i_plot_annotated_ec_after_soupx_rpca.out
#SBATCH --error=logs/09i_plot_annotated_ec_after_soupx_rpca.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs plots/09i_ec_after_soupx_rpca_annotated

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/09i_plot_annotated_ec_after_soupx_rpca.R
