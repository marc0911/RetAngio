#!/bin/bash -l
#SBATCH --job-name=retangio_12_ec_sct
#SBATCH --output=logs/12_prepare_sct_objects_ec_only.out
#SBATCH --error=logs/12_prepare_sct_objects_ec_only.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo
mkdir -p logs results/objects results/sct

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/12_prepare_sct_objects_ec_only.R
