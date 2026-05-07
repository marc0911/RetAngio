#!/bin/bash -l
#SBATCH --job-name=09m_tip_markers
#SBATCH --output=logs/09m_findmarkers_rpca_tipcell_local_clusters_%j.out
#SBATCH --error=logs/09m_findmarkers_rpca_tipcell_local_clusters_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo
mkdir -p logs

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/09m_findmarkers_rpca_tipcell_local_clusters.R
