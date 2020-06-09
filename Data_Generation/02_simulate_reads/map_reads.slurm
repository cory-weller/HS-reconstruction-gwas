#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 8G
#SBATCH -t 0-12:00:00
#SBATCH -p standard
#SBATCH --account berglandlab


module load bwa/0.7.17
module load samtools/1.9
module load gatk/4.0.0.0

population_stem=${1}

chromosome=${2}

first_ind=${3}

last_ind=${4}

generation=${SLURM_ARRAY_TASK_ID}

population=${population_stem}_F${generation}

projectDir=${PWD%02_simulate_reads}

tmpDir="/dev/shm/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

# tmpDir="/scratch/caw5cv/dev/shm/"

function mapIndividual() {
    local filestem=${1}
    local tmpDir=${2}
    mkdir -p ${tmpDir}/${filestem} && cd ${tmpDir}/${filestem}
    
    # run PEAR to assemble paired end reads
    ${projectDir}/bin/pear-0.9.11-linux-x86_64/bin/pear -k \
    -f ${projectDir}/02_simulate_reads/${population}/${filestem}.F.fq.bz2 \
    -r ${projectDir}/02_simulate_reads/${population}/${filestem}.R.fq.bz2 \
    -o ${filestem}
    
    # map assembled reads
    bwa mem ${projectDir}/input_data/dgrp2.reference.fasta \
    ${filestem}.assembled.fastq | \
    samtools view -b \
    -o ${filestem}.assembled.bam -
    
    # map unassembled reads
    bwa mem ${projectDir}/input_data/dgrp2.reference.fasta \
    ${filestem}.unassembled.forward.fastq ${filestem}.unassembled.reverse.fastq | \
    samtools view -b \
    -o ${filestem}.unassembled.bam -
    
    # merge bam files with samtools
    samtools merge ${filestem}.merged.bam ${filestem}.assembled.bam ${filestem}.unassembled.bam
    
    # Sort bam file
    samtools sort ${filestem}.merged.bam > ${filestem}.sorted.bam
    
    # Add read groupsand index .bam file
    echo "adding read groups"
    gatk AddOrReplaceReadGroups \
    --INPUT ${filestem}.sorted.bam \
    --OUTPUT ${filestem}.bam \
    --RGLB library1 \
    --RGPL illumina \
    --RGPU platform1 \
    --RGSM ${filestem}
    
    # index final output
    samtools index ${filestem}.bam
    
    # write final output to permenant disk location
    cp ${filestem}.bam* ${projectDir}/02_simulate_reads/${population}/
    
    # clean up
    cd ${projectDir}/02_simulate_reads && \
    rm -rf ${tmpDir}/${filestem}
}

export -f mapIndividual

for individual_n in $(seq ${first_ind} ${last_ind}); do
  filestem=${individual_n}.${chromosome}
  echo "mapping ${filestem}"
  mapIndividual ${filestem} ${tmpDir}
done

rm -rf ${tmpDir}
