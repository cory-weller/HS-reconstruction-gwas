#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=5
#SBATCH --mem 32G
#SBATCH -t 0-8:00:00
#SBATCH -p standard
#SBATCH --account berglandlab

module load gcc
module load R/3.5.1
module load tabix

count=${SLURM_ARRAY_TASK_ID}

population=${1}
replicates=${2}
n_loci=${3}
effect_size=${4}

folder_name="${population}_${count}_${n_loci}_${effect_size}"

mkdir -p ${folder_name} && cd ${folder_name}

echo "running liability.R ${population} ${replicates} ${n_loci} ${effect_size}"

Rscript ../liability.R ${population} ${replicates} ${n_loci} ${effect_size}

zip -Tm ../${folder_name}.zip * && cd ../ && rm -rf ${folder_name}
