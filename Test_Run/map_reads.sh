#!/usr/bin/env bash

REFGENOME_FILENAME="dmel-all-chromosome-r5.57.fasta"
SORTED_VCF="dgrp2.sorted.vcf.gz"

for FILENAME in *_R1_001.fastq; do
    FILESTEM=${FILENAME%_R1_001.fastq}
    ./map.sh ${FILESTEM} ${REFGENOME_FILENAME} ${SORTED_VCF}
done