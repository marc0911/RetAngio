#!/bin/bash -l
#SBATCH --job-name=10ya_sample_colors
#SBATCH --output=logs/10ya_fix_overall_sample_contribution_stacked_colorlocked.%j.out
#SBATCH --error=logs/10ya_fix_overall_sample_contribution_stacked_colorlocked.%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/10ya_fix_overall_sample_contribution_stacked_colorlocked.R
