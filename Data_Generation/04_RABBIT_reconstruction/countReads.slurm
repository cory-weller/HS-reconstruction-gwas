#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 8G
#SBATCH -t 0-12:00:00
#SBATCH -p standard
#SBATCH --account berglandlab

# submit as
# sbatch --array=1,2,5,50%1 countReaads.slurm

module load gparallel
module load gatk

# edit if different number of total individuals
N_individuals=1000

population="hybrid_swarm_32_F${SLURM_ARRAY_TASK_ID}"

# tmpWorkDir=/dev/shm//${USER}/${JOBID}/${population}
# tmpWorkDir="/scratch/${USER}/temp/"

projectDir="/scratch/${USER}/genome-reconstruction-revision/"
outputDir="/scratch/${USER}/genome-reconstruction-revision/04_RABBIT/${population}/readcounts/"
mkdir -p ${outputDir}

# harp parameters
# window_step=100000
# window_width=1000000
chromosome="2L"
# chrLength="23011544"
# harp="/scratch/caw5cv/genome-reconstruction-revision/04_RABBIT/harp"
referenceGenome="${projectDir}/input_data/dgrp2.reference.fasta"
bamFolder="${projectDir}/02_simulate_reads/${population}/"
variants_file="${projectDir}/input_data/variants_het_2L.vcf"
# priorsGzFile="${projectDir}/04_RABBIT/2L.priors.csv.gz"

#echo creating directory "${tmpWorkDir}"
#mkdir -p ${tmpWorkDir} && cd ${tmpWorkDir} && gunzip -c ${priorsGzFile} > priors.csv
#priors="${PWD}/priors.csv"

# for each individual in population...

function getReadCounts {
  ind=${1}
  bamFolder=${2}
  chromosome=${3}
  #tmpWorkDir=${4}
  outputDir=${4}
  #chrLength=${6}
  referenceGenome=${5}
  #priors=${8}
  #window_step=${9}
  #window_width=${10}
  #harp=${11}
  variants_file=${6}

  bamPath="${bamFolder}/${ind}.${chromosome}.bam"

  # make tmp work directory
  #echo creating directory "${tmpWorkDir}/${ind}"
  #mkdir -p "${tmpWorkDir}/${ind}" && cd "${tmpWorkDir}/${ind}"

  echo "calculating readcounts for ${ind}"
  
  gatk --java-options '-Xmx6G' ASEReadCounter \
  --reference ${referenceGenome} \
  --input ${bamPath} \
  --output ${outputDir}/$ind.${chromosome}.readcounts \
  --variant ${variants_file}
}

export -f getReadCounts

# single run:
parallel -j 1 getReadCounts ::: $(seq 1 $N_individuals) ::: ${bamFolder} ::: ${chromosome} ::: ${outputDir} ::: ${referenceGenome} ::: ${variants_file}
