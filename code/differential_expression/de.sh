#!/bin/bash
#SBATCH --job-name=Atlas_DE
#SBATCH --output=DE_%j.out
#SBATCH --mem=200G
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00

# Load both R and the Bioconductor bundle
module load r/4.5.0
module load r-bundle-bioconductor/3.21

Rscript de.R
