#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=5
#SBATCH --mem 8G
#SBATCH -t 0-04:00:00
#SBATCH -p standard
#SBATCH --account berglandlab


module load gcc
module load R/3.5.1

Rscript roc.R
