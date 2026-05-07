#!/bin/bash -l
#SBATCH --job-name=retangio_08_unint
#SBATCH --output=logs/08_unint.out
#SBATCH --error=logs/08_unint.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/08_unintegrated_clustering_umap.R
