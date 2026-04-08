#!/bin/bash
#SBATCH --job-name=Elbow_Check
#SBATCH --output=elbow_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=02:00:00

module load r/4.5.0
module load r-bundle-bioconductor/3.21

Rscript check_dimensions.R
