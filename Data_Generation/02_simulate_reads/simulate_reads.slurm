#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 8G
#SBATCH -t 0-12:00:00
#SBATCH -p standard
#SBATCH --account berglandlab


population_stem=${1}

chromosome=${2}

first_ind=${3}

last_ind=${4}

generation=${SLURM_ARRAY_TASK_ID}

population=${population_stem}_F${generation}

python simulate_reads.py ${population} ${chromosome} ${first_ind} ${last_ind}
