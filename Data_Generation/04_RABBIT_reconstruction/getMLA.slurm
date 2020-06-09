#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem 32G
#SBATCH -t 0-12:00:00
#SBATCH -p standard
#SBATCH --account berglandlab

module load gcc
module load R/3.5.1
module load gparallel

population_stem="${1}"
chromosome="${2}"
first_ind="${3}"
last_ind="${4}"

generation="${SLURM_ARRAY_TASK_ID}"

population="${population_stem}_F${generation}"

mkdir -p "${population}/MLA/"

parallel -j 4 Rscript ./getMLA.R ::: ${population} ::: $(seq ${first_ind} ${last_ind}) ::: ${chromosome}


# for ind_n in $(seq 1 ${N_individuals}); do
#   Rscript ./getMLA.R ${population} ${ind_n} ${chromosome}
# done
