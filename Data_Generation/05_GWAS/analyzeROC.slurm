#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --mem 32G
#SBATCH -t 0-2:00:00
#SBATCH -p standard
#SBATCH --account berglandlab

module load gcc
module load R/3.5.1

Rscript analyzeROC.R
