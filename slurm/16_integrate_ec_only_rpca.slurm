#!/bin/bash -l
#SBATCH --job-name=retangio_16_ec_rpca
#SBATCH --output=logs/16_integrate_ec_only_rpca.out
#SBATCH --error=logs/16_integrate_ec_only_rpca.err
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo
mkdir -p logs results/objects results/integration

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/16_integrate_ec_only_rpca.R
