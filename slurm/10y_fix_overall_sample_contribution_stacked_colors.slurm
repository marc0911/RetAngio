#!/bin/bash
#SBATCH --job-name=10y_overall_stack_fix
#SBATCH --output=logs/10y_fix_overall_sample_contribution_stacked_colors.%j.out
#SBATCH --error=logs/10y_fix_overall_sample_contribution_stacked_colors.%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

module purge
module load r-light

export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript scripts/10y_fix_overall_sample_contribution_stacked_colors.R
