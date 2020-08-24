#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem 24G
#SBATCH -t 0-20:00:00
#SBATCH -p standard
#SBATCH --account berglandlab

# submit as
# sbatch --array=1,2,5,50%1 doHarp.slurm

module load gparallel

# edit if different number of total individuals
N_individuals=1000

population="hybrid_swarm_32_F${SLURM_ARRAY_TASK_ID}"

tmpWorkDir=/dev/shm//${USER}/${JOBID}/${population}
# tmpWorkDir="/scratch/${USER}/temp/"

projectDir="/scratch/${USER}/genome-reconstruction-revision/"
outputDir="/scratch/${USER}/genome-reconstruction-revision/04_RABBIT/${population}/HARP_likelihood/"
mkdir -p ${outputDir}

declare -A lengths='([3L]="24543557" [X]="22422827" [3R]="27905053" [2R]="21146708" [2L]="23011544" )'

# harp parameters
window_step=100000
window_width=100000
chromosome="2L"
chrLength="23011544"
harp="/scratch/caw5cv/genome-reconstruction-revision/04_RABBIT/harp"
referenceGenome="${projectDir}/input_data/dgrp2.reference.fasta"
bamFolder="${projectDir}/02_simulate_reads/${population}/"
priorsGzFile="${projectDir}/04_RABBIT/${population}/${chromosome}.subset.priors.csv.gz"

if [ ! -f "${priorsGzFile}" ]; then
  Rscript ./subsetPriors.R ${population} ${chromosome} && \
  sed -i 's/$/,/g' ${population}/${chromosome}.subset.priors.csv && \
  gzip ${population}/${chromosome}.subset.priors.csv
fi

echo creating directory "${tmpWorkDir}"
mkdir -p ${tmpWorkDir} && cd ${tmpWorkDir} && gunzip -c ${priorsGzFile} > priors.csv
priors="${PWD}/priors.csv"

# for each individual in population...

function getHarpFreqs {
  ind=${1}
  bamFolder=${2}
  chromosome=${3}
  tmpWorkDir=${4}
  outputDir=${5}
  chrLength=${6}
  referenceGenome=${7}
  priors=${8}
  window_step=${9}
  window_width=${10}
  harp=${11}

  bamPath="${bamFolder}/${ind}.${chromosome}.bam"

  # make tmp work directory
  echo creating directory "${tmpWorkDir}/${ind}"
  mkdir -p "${tmpWorkDir}/${ind}" && cd "${tmpWorkDir}/${ind}"

  # Run harp like
  echo running harp like
  $harp like \
  --bam $bamPath \
  --region $chromosome:1-$chrLength \
  --refseq $referenceGenome \
  --snps $priors \
  --stem $ind.$chromosome

  # Run harp freq
  echo running harp freq
  $harp freq \
  --bam $bamPath \
  --region $chromosome:1-$chrLength \
  --refseq $referenceGenome \
  --snps $priors \
  --stem $ind.$chromosome \
  --window_step $window_step \
  --window_width $window_width \
  --em_min_freq_cutoff 0.0001

  # Cleanup and move out
  echo cleaning up
  cp "${ind}.${chromosome}.freqs" ${outputDir} && \
  cd ${tmpWorkDir} && \
  rm -rf ${tmpWorkDir}/${ind}
  echo done
}

export -f getHarpFreqs

# single run:
parallel -j 4 getHarpFreqs ::: $(seq 1 $N_individuals) ::: $bamFolder ::: $chromosome ::: $tmpWorkDir ::: $outputDir ::: $chrLength ::: $referenceGenome ::: $priors ::: $window_step ::: $window_width ::: $harp

rm -rf ${tmpWorkDir}
