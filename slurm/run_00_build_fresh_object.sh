#!/bin/bash -l
#SBATCH --job-name=retangio_00_build
#SBATCH --output=logs/00_build.out
#SBATCH --error=logs/00_build.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# Conservative defaults for stage 0.
# Adjust time / CPU / memory after testing on your real data.

set -euo pipefail

mkdir -p logs

module load R

Rscript R/00_build_fresh_object.R config/config.yml
