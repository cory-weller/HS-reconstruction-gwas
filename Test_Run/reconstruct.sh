#!/usr/bin/env bash

BAM_FILESTEM=${1%.final.bam}
BAM_FILENAME=$(realpath ${1})
REFGENOME_FILENAME=$(realpath ${2})
UTILS_SIF=$(realpath ../utils.sif)
HARP_SIF=$(realpath ../harp.sif)
TMP="/home/cory/tmpStorage"
ORIG_DIR=$(pwd)

# cycles through every chromosome within a .bam file to generate .freqs output

# Read chromosome lengths into associative array from the fasta index
declare -A LENGTHS
while read -r CHROM LENGTH; do
    LENGTHS[${CHROM}]=${LENGTH}
done < <(awk '{print $1,$2}' ${REFGENOME_FILENAME}.fai)

function getHarpFreqs {
    HARP_STEP=100000
    HARP_WIDTH=100000
    BAM_FILENAME=${1}
    BAM_FILESTEM=${2}
    CHROMOSOME=${3}
    CHROMOSOME_LENGTH=${4}
    REFGENOME_FILENAME=${5}
    ORIG_DIR=${6}
    TMP=${7}
    TMP_WORK_DIR="${TMP}/${BAM_FILESTEM}/${CHROMOSOME}"

    HARP_CSV_FILENAME="${ORIG_DIR}/${CHROMOSOME}.harp.csv"

    # make tmp work directory
    echo creating directory "${TMP_WORK_DIR}"
    mkdir -p "${TMP_WORK_DIR}" && cd "${TMP_WORK_DIR}/"

    # Clean up tmp directory on unexpected exit
    trap 'rm -rf "${TMP_WORK_DIR}"' EXIT

    # Run harp like
    echo running harp like
        singularity exec ${HARP_SIF} harp like \
        --bam ${BAM_FILENAME} \
        --region ${CHROMOSOME}:1-${CHROMOSOME_LENGTH} \
        --refseq ${REFGENOME_FILENAME} \
        --snps ${HARP_CSV_FILENAME} \
        --stem ${BAM_FILESTEM}.${CHROMOSOME}

    # Run harp freq
    echo running harp freq
    singularity exec ${HARP_SIF} harp freq \
        --bam ${BAM_FILENAME} \
        --region ${CHROMOSOME}:1-${CHROMOSOME_LENGTH} \
        --refseq ${REFGENOME_FILENAME} \
        --snps ${HARP_CSV_FILENAME} \
        --stem ${BAM_FILESTEM}.${CHROMOSOME} \
        --window_step ${HARP_STEP} \
        --window_width ${HARP_WIDTH} \
        --em_min_freq_cutoff 0.0001

    # Cleanup and move out
    echo cleaning up
    echo $(ls)
    cp "${BAM_FILESTEM}.${CHROMOSOME}.freqs" ${ORIG_DIR} && \
    cd ${ORIG_DIR} && \
    rm -rf ${TMP_WORK_DIR}
    echo done
}

export -f getHarpFreqs

# Iterate through chromosomes with extant *.harp.csv file
for FILENAME in *.harp.csv; do
    CHROMOSOME=${FILENAME%.harp.csv}
    CHROMOSOME_LENGTH=${LENGTHS[${CHROMOSOME}]}

    # Get HARP freqs
    if [[ ! -f "${BAM_FILESTEM}.${CHROMOSOME}.freqs" ]]; then
        echo "Running HARP on ${CHROMOSOME} for ${BAM_FILESTEM}"
        getHarpFreqs "${BAM_FILENAME}" "${BAM_FILESTEM}" "${CHROMOSOME}" "${CHROMOSOME_LENGTH}" "${REFGENOME_FILENAME}" "${ORIG_DIR}" "${TMP}"
    else
        echo "${BAM_FILESTEM}.${CHROMOSOME}.freqs already exists"
    fi

    # Get MLA (write .mla file) for this chromosome
    if [[ ! -f "${BAM_FILESTEM}.${CHROMOSOME}.mla" ]]; then
        echo "Calculating most likely ancestors for ${CHROMOSOME} for ${BAM_FILESTEM}"
        singularity exec ${UTILS_SIF} Rscript get_MLA.R ${BAM_FILESTEM} ${CHROMOSOME}
    else
        echo "${BAM_FILESTEM}.${CHROMOSOME}.mla already exists"
    fi

done


