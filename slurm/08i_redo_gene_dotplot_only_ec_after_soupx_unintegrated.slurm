#!/bin/bash -l
#SBATCH --job-name=retangio_08i_dot
#SBATCH --output=logs/08i_redo_gene_dotplot_only_ec_after_soupx_unintegrated.out
#SBATCH --error=logs/08i_redo_gene_dotplot_only_ec_after_soupx_unintegrated.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G

set -euo pipefail

cd /work/PRTNR/CHUV/HOJG/mschwab2/retangio/repo

mkdir -p logs plots/08i_ec_after_soupx_unintegrated_signatures_and_genes

module purge
module load r-light
export R_LIBS_USER=/work/PRTNR/CHUV/HOJG/mschwab2/retangio/Rlibs

Rscript R/08i_redo_gene_dotplot_only_ec_after_soupx_unintegrated.R
