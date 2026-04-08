#!/bin/bash
#SBATCH --job-name=Atlas_UMAP_Final
#SBATCH --output=UMAP_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=350G
#SBATCH --time=05:00:00

module load r/4.5.0
module load r-bundle-bioconductor/3.21

Rscript UMAP.R
