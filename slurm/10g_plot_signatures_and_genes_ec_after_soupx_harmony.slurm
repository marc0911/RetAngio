#!/bin/bash -l
#SBATCH --job-name=retangio_10g_plots
#SBATCH --output=logs/10g_plot_signatures_and_genes_ec_after_soupx_harmony.out
#SBATCH --error=logs/10g_plot_signatures_and_genes_ec_after_soupx_harmony.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs plots/10g_ec_after_soupx_harmony_signatures_and_genes

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/10g_plot_signatures_and_genes_ec_after_soupx_harmony.R
