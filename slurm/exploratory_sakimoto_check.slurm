#!/bin/bash -l
#SBATCH --job-name=retangio_sakimoto_check
#SBATCH --output=logs/exploratory_sakimoto_check.out
#SBATCH --error=logs/exploratory_sakimoto_check.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo
mkdir -p logs results/exploratory_sakimoto plots/exploratory_sakimoto

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/exploratory_sakimoto_check.R
