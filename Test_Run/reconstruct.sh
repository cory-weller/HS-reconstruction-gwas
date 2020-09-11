#!/usr/bin/env bash

BAM_FILESTEM=${1%.final.bam}
BAM_FILENAME="$(pwd)/${1}"
REFGENOME_FILENAME="$(pwd)/${2}"
N_GENERATIONS=${3}
UTILS_SIF="$(pwd)/../utils.sif"
HARP_SIF="$(pwd)/../harp.sif"
TMP="/dev/shm/caw5cv/tmp/"
ORIG_DIR=$(pwd)

RABBIT_MAX_FOUNDERS="16"
RABBIT_MAX_SITES="5000"



# cycles through every chromosome within a .bam file to generate .freqs output

# Read chromosome lengths into associative array from lenghts.txt
declare -A LENGTHS
while read -r CHROM LENGTH; do
    LENGTHS[${CHROM}]=${LENGTH}
done < lengths.txt

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
for CHROMOSOME in ${!LENGTHS[@]}; do
     cd ${ORIG_DIR}
    CHROMOSOME_LENGTH=${LENGTHS[${CHROMOSOME}]}
    HAPLOTYPES_FILE="${BAM_FILESTEM}.${CHROMSOME}.haplotypes"
    RABBIT_CSV="${BAM_FILESTEM}.${CHROMSOME}.RABBIT.out.csv"

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

    # Count Reads
    if [[ ! -f ${BAM_FILESTEM}.${CHROMOSOME}.readcounts ]]; then
    echo "Running ASEReadCounter on ${CHROMOSOME} for ${BAM_FILESTEM}"
    singularity exec ${UTILS_SIF} java -jar /opt/gatk4.jar ASEReadCounter \
        --reference ${REFGENOME_FILENAME} \
        --input ${BAM_FILENAME} \
        --output ${BAM_FILESTEM}.${CHROMOSOME}.readcounts \
        --variant ${CHROMOSOME}.variants.het.vcf
    else
        echo "${BAM_FILESTEM}.${CHROMOSOME}.readcounts already exists"
    fi

    # Build RABBIT input
    if [[ ! -f ${BAM_FILESTEM}.${CHROMOSOME}.RABBIT.in ]]; then
        echo "Generating ${BAM_FILESTEM}.${CHROMOSOME}.RABBIT.in"
        singularity exec ${UTILS_SIF} Rscript build_RABBIT_input.R ${BAM_FILESTEM} ${CHROMOSOME} ${RABBIT_MAX_FOUNDERS} ${RABBIT_MAX_SITES}
    else
        echo "${BAM_FILESTEM}.${CHROMOSOME}.RABBIT.in already exists"
    fi

    # If F1 population, increment generation by 1 to allow heterozygous RABBIT imputation
    if [[ ${N_GENERATIONS} == "1" ]]; then
        let RABBIT_GENERATIONS=${N_GENERATIONS}+1
    else
        let RABBIT_GENERATIONS=${N_GENERATIONS}
    fi

    if [[ ! -f ${BAM_FILESTEM}.${CHROMOSOME}.RABBIT.m ]]; then
        echo "Generating ${BAM_FILESTEM}.${CHROMOSOME}.RABBIT.m"
        # Generate mathematica script for RABBIT
        CURRENT_DIR=${ORIG_DIR}
        RABBIT_DIR="${ORIG_DIR}/RABBIT/RABBIT_Packages/"
        RABBIT_eps="0.005"
        RABBIT_epsF="0.005"
        RABBIT_MODEL="jointModel"
        RABBIT_EST_FUN="origViterbiDecoding"

        singularity exec ${UTILS_SIF} python3 - <<EOF > ${BAM_FILESTEM}.${CHROMOSOME}.RABBIT.m
print("""SetDirectory["%s"]""" % "${RABBIT_DIR}")
print("""Needs["MagicReconstruct\`"]""")
print("""SetDirectory["%s"]""" % "${CURRENT_DIR}/")
print("""popScheme = Table["RM1-E", {%s}]""" % "${RABBIT_GENERATIONS}")
print('epsF = %s' % "${RABBIT_epsF}")
print('eps = %s' % "${RABBIT_eps}")
print('model = "%s"' % "${RABBIT_MODEL}")
print('estfun = "%s"' % "${RABBIT_EST_FUN}")
print('inputfile = "%s"' % "${BAM_FILESTEM}.${chromosome}.RABBIT.in")
print('resultFile = "%s.txt"' % "${BAM_FILESTEM}.${chromosome}.RABBIT.out")
print("""magicReconstruct[inputfile, model, RABBIT_epsF, eps, popScheme, resultFile, HMMMethod -> estfun, PrintTimeElapsed -> True]""")
print('summaryFile = StringDrop[resultFile, -4] <> ".csv"')
print('saveAsSummaryMR[resultFile, summaryFile]')
print('Exit')
EOF

    else
        echo "${BAM_FILESTEM}.${CHROMOSOME}.RABBIT.m already exists"
    fi
done



if [ -f "${RABBIT_CSV}" ]; then
    echo "${RABBIT_CSV}" already exists
    echo skipping RABBIT for this individual
    continue
else
    echo "running mathematica"
    math -noprompt -script "${BAM_FILESTEM}.${CHROMOSOME}.RABBIT.m"
fi


# Print haplotypes file header
echo -e "chromosome\tstart\tstop\tpar1\tpar2" > "${HAPLOTYPES_FILE}"

# if nonrecombinant (i.e. no rabbit-generated output, only manual output),
# parse manual output to make nonrecombinant haplotype map file
N_LINES=$(wc -l ${RABBIT_CSV} | awk '{print $1}')
if [ "${N_LINES}" == 1 ]; then
    ID=$(awk '{print $3}' ${RABBIT_CSV})
    echo -e "${CHROMOSOME}\t1\t25000000\t${ID}\t${ID}" >> "${HAPLOTYPES_FILE}"
else
    singularity exec ${UTILS_SIF} python3 convert_to_haplotypes.py ${RABBIT_CSV} >> "${HAPLOTYPES_FILE}"
fi

done
