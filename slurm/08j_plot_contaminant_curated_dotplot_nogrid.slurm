#!/bin/bash -l
#SBATCH --job-name=08j_contam_dot
#SBATCH --output=logs/08j_plot_contaminant_curated_dotplot_nogrid.%j.out
#SBATCH --error=logs/08j_plot_contaminant_curated_dotplot_nogrid.%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/08j_plot_contaminant_curated_dotplot_nogrid.R
