#!/bin/bash -l
#SBATCH --job-name=retangio_09f_sign
#SBATCH --output=logs/09f_signature_scoring_ec_after_soupx_rpca.out
#SBATCH --error=logs/09f_signature_scoring_ec_after_soupx_rpca.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs results/signatures results/objects

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/09f_signature_scoring_ec_after_soupx_rpca.R
