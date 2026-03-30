#!/usr/bin/env bash
#SBATCH --job-name=retangio_00_build
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

# Conservative defaults above are usually safe for stage 0,
# but you may need to adjust CPU/memory/time for your cluster and dataset size.

mkdir -p logs

module load R

Rscript R/00_build_fresh_object.R config/config.yml
