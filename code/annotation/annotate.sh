#!/bin/bash
#SBATCH --job-name=Atlas_Annotate
#SBATCH --output=annotate_%j.out
#SBATCH --mem=500G
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00

# Load both R and the Bioconductor bundle
module load r/4.5.0
module load r-bundle-bioconductor/3.21

Rscript annotate.R
