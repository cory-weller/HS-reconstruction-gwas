#!/usr/bin/env bash

if ! command -v realpath &> /dev/null; then
    echo "command realpath could not be found. Install coreutils?"
    exit
fi

source PARAMETERS.config


if [ "${REFGENOME_URL: -3}" == ".gz" ]; then
    REFGENOME_ZIPPED=$(basename "${REFGENOME_URL}")
    REFGENOME="${REFGENOME_ZIPPED%.gz}"
elif [ "${REFGENOME_URL: -6}" == ".fasta" ] || [ "${REFGENOME_URL: -3}" == ".fa" ]; then
    REFGENOME=$(basename "${REFGENOME_URL}")
fi

if [ "${VCF_URL: -3}" == ".gz" ]; then
    VCF_ZIPPED=$(basename "${VCF_URL}")
    VCF="${VCF_ZIPPED%.gz}"
elif [ "${VCF_URL: -4}" == ".vcf" ]; then
    VCF=$(basename "${VCF_URL}")
    VCF_ZIPPED="${VCF}.gz"
fi


# Create reference sequence dictionary
if [ ! -f "${REFGENOME%.fasta}".dict ]; then
    echo "Building reference sequence dictionary"
    singularity exec ${UTILS_SIF} java -jar /opt/gatk4.jar CreateSequenceDictionary \
        -R "${REFGENOME}" \
        -O "${REFGENOME%.fasta}".dict
else
    echo "Reference sequence dictionary ${REFGENOME}.fasta.dict already generated..."
fi

# Prepare SAMTOOLS index
if [ ! -f "${REFGENOME}.fai" ]; then
    echo "Building samtools index"
    singularity exec ${HARP_SIF} samtools faidx ${REFGENOME}
else
    echo "Samtools index ${REFGENOME}.fai already generated..."
fi

# Prepare BWA index
if [ ! -f "${REFGENOME}.amb" ]; then
    echo "Building BWA index"
    singularity exec ${UTILS_SIF} bwa index ${REFGENOME}
else
    echo "BWA index already generated..."
fi

# Generate lengths.txt file containing chromosome lengths
declare -A LENGTHS
if [ -f lengths.txt ]; then rm lenghts.txt; fi
while read -r CHROM LENGTH; do
    if printf '%s\n' "${CHROMOSOMES[@]}" | grep -q -P "^${CHROM}$"; then
        echo ${CHROM} ${LENGTH} >> lengths.txt
    fi
done < <(awk '{print $1,$2}' ${REFGENOME}.fai)

# If incorrect md5 sum; remove corrupt file
if [ -f ${VCF_ZIPPED} ] && [ ! $(md5sum ${VCF_ZIPPED} | awk '{print $1}') == ${VCF_BGZIP_MD5} ]; then
    echo "VCF has incorrrect md5 sum. Re-download. Exiting."
    exit 1
fi

# Sort VCF
if [[ ! -f ${VCF_ZIPPED%.vcf.gz}.sorted.vcf ]]; then
    echo "Sorting VCF..."
    # dgrp2.sorted.vcf md5sum: dac0ddb8da32d9e4b68458b7aff28caf
    singularity exec ${UTILS_SIF} java -Xmx2G -jar /opt/gatk4.jar SortVcf \
            -I ${VCF_ZIPPED} \
            -O ${VCF_ZIPPED%.vcf.gz}.sorted.vcf \
            --SEQUENCE_DICTIONARY ${REFGENOME%.fasta}.dict \
            --MAX_RECORDS_IN_RAM 100000
    echo "BGZIP on sorted VCF..."
    singularity exec ${UTILS_SIF} bgzip -f ${VCF_ZIPPED%.vcf.gz}.sorted.vcf
    # dgrp2.sorted.vcf.gz md5sum: e8364935cb63d5fc945a072321d8631c
fi

# Prepare VCF subset without indels
singularity exec ${UTILS_SIF} awk -v file=${VCF_ZIPPED%.vcf.gz}.sorted.noIndel.vcf -F "\t" '{
if($0 ~ /^#/ || $3 ~ /_SNP/)
        print >> file}' < <(zcat ${VCF_ZIPPED%.vcf.gz}.sorted.vcf.gz)


# Prepare VCF subset with only indels
singularity exec ${UTILS_SIF} awk -v file=${VCF_ZIPPED%.vcf.gz}.sorted.indel.vcf -F "\t" '{
if($0 ~ /^#/ || $3 ~ /_DEL/ || $3 ~ /_INS/)
        print >> file}' < <(zcat ${VCF_ZIPPED%.vcf.gz}.sorted.vcf.gz)


# dgrp2.sorted.indel.vcf md5: 0fe9296e7c1681697492c9175a43986f
# dgrp2.sorted.noIndel.vcf md5: c4793c71cf076964441944845ab8b921


# Format as samtools ranges, e.g. CHROM:START-STOP
# repeatmasker is 1-based, closed [start, end] 
# bed is    0-based, half-open [start-1, end)
# Format as 1-based .intervals and concatenate
for CHR in ${CHROMOSOMES[@]}; do
    cat <(awk 'BEGIN{OFS="\t"} NR > 3 {printf ("%s:%s-%s\n",$5,$6,$7)}' tmp_repetitive/${CHR}/chr${CHR}.fa.out) \
            <(awk 'BEGIN{OFS="\t"} {printf ("%s:%s-%s\n", $1,$2+1,$3)}' tmp_repetitive/trfMaskChrom/chr${CHR}.bed) | \
            sort -nk 2,3 | \
            sed "s/chr//g" >> repetitive.list
done && \
rm -rf tmp_repetitive/

# Index sorted vcf
singularity exec ${UTILS_SIF} java -jar /opt/gatk4.jar IndexFeatureFile \
    --input ${VCF%.vcf}.sorted.noIndel.vcf

# Exclude repeptitive intervals
singularity exec ${UTILS_SIF} java -Xmx4G -jar /opt/gatk4.jar SelectVariants \
    --variant    ${VCF%.vcf}.sorted.noIndel.vcf \
    --exclude-intervals repetitive.list \
    --output ${VCF%.vcf}.sorted.noIndel.noRep.vcf

# Generate HARP genotype csv files
singularity exec ${UTILS_SIF} python3 generate_harp_csv.py ${VCF%.vcf}.sorted.noIndel.noRep.vcf

# Generate HARP input haplotypes
for CHROMOSOME in ${CHROMOSOMES[@]}; do
    singularity exec ${UTILS_SIF} python3 generate_RABBIT_csv.py ${VCF%.vcf}.sorted.noIndel.noRep.vcf
done


# prepare heterozygous VCF file for ASEReadCounter
for CHROMOSOME in ${CHROMOSOMES[@]}; do
    awk -F "\t" -v CHR=${CHROMOSOME}    '
    BEGIN{OFS="\t";}
    {
    if ($1 ~ /^##/)
                    print $0;
    else if ($1 ~ /^#CHROM/)
                    print $1,$2,$3,$4,$5,$6,$7,$8,$9,"HET";
    else if ($1 ~ CHR)
                    print $1,$2,$3,$4,$5,$6,$7,$8,$9,"0/1";
    }' <(cat ${VCF%.vcf}.sorted.noIndel.noRep.vcf) > ${CHROMOSOME}.variants.het.vcf

    singularity exec ${UTILS_SIF} java -jar /opt/gatk4.jar IndexFeatureFile --input ${CHROMOSOME}.variants.het.vcf
done



exit 0