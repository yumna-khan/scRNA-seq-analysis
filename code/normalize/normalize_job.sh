#!/bin/bash
#SBATCH --job-name=IAV_Atlas_Proc
#SBATCH --output=logs_%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=12:00:00

module load r/4.5.0
module load r-bundle-bioconductor/3.21

Rscript run_normalization.R
