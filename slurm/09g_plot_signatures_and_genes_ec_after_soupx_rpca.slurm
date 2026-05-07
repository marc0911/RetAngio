#!/bin/bash -l
#SBATCH --job-name=retangio_09g_plots
#SBATCH --output=logs/09g_plot_signatures_and_genes_ec_after_soupx_rpca.out
#SBATCH --error=logs/09g_plot_signatures_and_genes_ec_after_soupx_rpca.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs plots/09g_ec_after_soupx_rpca_signatures_and_genes

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/09g_plot_signatures_and_genes_ec_after_soupx_rpca.R
