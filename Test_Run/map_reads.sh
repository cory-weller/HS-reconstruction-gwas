#!/usr/bin/env bash
# Run as: bash ./map_reads.sh ${FILE_ID} ${REFGENOME_FILENAME} ${SORTED_VCF}
# e.g.    bash ./map_reads.sh 4CB-1-82_S82 dmel-all-chromosome-r5.57.fasta dgrp2.sorted.vcf

UTILS_SIF="../utils.sif"
HARP_SIF="../harp.sif"
java8="/usr/lib/jvm/java-8-openjdk-amd64/bin/java" 
# java11="/usr/lib/jvm/java-11-openjdk-amd64/bin/java" Equivalent to `java`

FILE_ID=${1}
REFGENOME_FILENAME=${2}
SORTED_VCF=${3}

VCF_FILESTEM=${SORTED_VCF%.sorted.vcf}

# run PEAR
singularity exec ${UTILS_SIF} /opt/pear \
    -f "$FILE_ID"_R1_001.fastq \
    -r "$FILE_ID"_R2_001.fastq \
    -o "$FILE_ID"

# Map unassembled reads to reference genome
singularity exec ${UTILS_SIF} bwa mem \
    ${REFGENOME_FILENAME} \
    ${FILE_ID}.unassembled.forward.fastq ${FILE_ID}.unassembled.reverse.fastq | \
    singularity exec ${HARP_SIF} samtools view -Shu - | \
    singularity exec ${HARP_SIF} samtools sort - "$FILE_ID".unassembled.sorted

# Map assembled reads to reference genome
singularity exec ${UTILS_SIF} bwa mem \
    ${REFGENOME_FILENAME} \
    ${FILE_ID}.assembled.fastq | \
    singularity exec ${HARP_SIF} samtools view -Shu - | \
    singularity exec ${HARP_SIF} samtools sort - "$FILE_ID".assembled.sorted

# Merge assembled and unassembled bam files
singularity exec ${HARP_SIF} samtools merge \
    ${FILE_ID}.merged.bam \
    ${FILE_ID}.assembled.sorted.bam \
    ${FILE_ID}.unassembled.sorted.bam

# Sort bam file by coordinates
singularity exec ${UTILS_SIF} java -Xmx2G -jar /opt/gatk4.jar SortSam \
      --INPUT ${FILE_ID}.merged.bam \
      --OUTPUT ${FILE_ID}.merged.sorted.bam \
      --SORT_ORDER coordinate

# Create bam file index
singularity exec ${HARP_SIF} samtools index ${FILE_ID}.merged.sorted.bam

# Mark duplicates with picard
singularity exec ${UTILS_SIF} java -Xmx2G -jar /opt/gatk4.jar MarkDuplicates \
    --INPUT ${FILE_ID}.merged.sorted.bam \
    --OUTPUT ${FILE_ID}.merged.markdup.bam \
    --METRICS_FILE ${FILE_ID}.duplicate_metrics.txt \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500

# Add required read groups
singularity exec ${UTILS_SIF} java -Xmx2G -jar /opt/gatk4.jar AddOrReplaceReadGroups \
    --INPUT ${FILE_ID}.merged.markdup.bam \
    --OUTPUT ${FILE_ID}.merged.markdup.rg.bam \
    --RGID ${FILE_ID} \
    --RGLB lib1 \
    --RGPL illumina \
    --RGPU unit1 \
    --RGSM ${FILE_ID} \
    --SORT_ORDER coordinate \
    --CREATE_INDEX true

# Create indel interval target file
singularity exec ${UTILS_SIF} ${java8} -Xmx2G -jar /opt/gatk3.jar \
    -T RealignerTargetCreator \
    -R ${REFGENOME_FILENAME} \
    -I ${FILE_ID}.merged.markdup.rg.bam \
    -o ${FILE_ID}.intervals

# Reorder BAM
singularity exec ${UTILS_SIF} java -Xmx2G -jar /opt/gatk4.jar ReorderSam \
    -I ${FILE_ID}.merged.markdup.rg.bam \
    -O ${FILE_ID}.merged.markdup.rg.reordered.bam \
    -R ${REFGENOME_FILENAME} \
    --SEQUENCE_DICTIONARY ${REFGENOME_FILENAME%.fasta}.dict \
    -CREATE_INDEX true

# Run Indel Realignment
singularity exec ${UTILS_SIF} ${java8} -jar /opt/gatk3.jar \
    -T IndelRealigner \
    -R ${REFGENOME_FILENAME} \
    -I ${FILE_ID}.merged.markdup.rg.reordered.bam \
    -known ${VCF_FILESTEM}.sorted.indel.vcf \
    -targetIntervals ${FILE_ID}.intervals \
    -o ${FILE_ID}.merged.markdup.rg.realign.bam

# Clean Up
