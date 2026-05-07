#!/bin/bash -l
#SBATCH --job-name=retangio_09h_annot
#SBATCH --output=logs/09h_annotate_ec_after_soupx_rpca.out
#SBATCH --error=logs/09h_annotate_ec_after_soupx_rpca.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs results/annotation results/objects

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/09h_annotate_ec_after_soupx_rpca.R
