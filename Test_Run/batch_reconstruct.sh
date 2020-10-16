#!/usr/bin/env bash

module purge

source ./PARAMETERS.config

ls *.final.bam > bamlist.txt

N_FILES=$(wc -l bamlist.txt | awk '{print $1}')

sbatch ${SLURM_OPTS} --array=1-"${N_FILES}"%3 ./reconstruct.sh
