#!/bin/bash -l
#SBATCH --job-name=09j_rpca_abundance
#SBATCH --output=logs/09j_plot_annotation_abundance_ec_after_soupx_rpca_%j.out
#SBATCH --error=logs/09j_plot_annotation_abundance_ec_after_soupx_rpca_%j.err
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/09j_plot_annotation_abundance_ec_after_soupx_rpca.R
