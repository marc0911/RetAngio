#!/usr/bin/env bash
#SBATCH --job-name=retangio_00_build
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

# Conservative defaults for stage 0.
# Adjust resources based on data size and cluster policy.

mkdir -p logs

module load R

Rscript R/00_build_fresh_object.R config/config.yml
