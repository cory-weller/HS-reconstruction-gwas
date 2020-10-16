
if ! command -v realpath &> /dev/null; then
    echo "command realpath could not be found. Install coreutils?"
    exit
fi

source PARAMETERS.config

if [ ! -f "${HARP_SIF}" ]; then
    singularity pull --name ${HARP_SIF} shub://cory-weller/HS-reconstruction-gwas:harp
fi

if [ ! -f "${UTILS_SIF}" ]; then
    singularity pull --name ${UTILS_SIF} shub://cory-weller/HS-reconstruction-gwas:utils
fi

if [[ ! -f ${RECOMBINATION_MAP_FILENAME} ]]; then
    wget -O ${RECOMBINATION_MAP_FILENAME} ${RECOMBINATION_MAP_URL}
fi

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



# download bam files
# wget -O 73.test_run.bam "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118194&authkey=AJNA5CF5UX1bsLI"
# wget -O 77.test_run.bam "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118191&authkey=AENQGxjwtHH3BD0"
# wget -O 78.test_run.bam "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118190&authkey=AHNAMKodHtgt0e0"
# wget -O 79.test_run.bam "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118192&authkey=ABxLTZCWwwuyjCY"
# wget -O 81.test_run.bam "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118195&authkey=AFdak4Jwyei0NzE"
# wget -O 82.test_run.bam "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118193&authkey=APj6cjxPbLdgl6c"

# Download low coverage reads
wget -O low_coverage_reads.zip "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118197&authkey=AAJg8Dw7UluE1w8" && \
singularity exec ${UTILS_SIF} unzip low_coverage_reads.zip && \
rm low_coverage_reads.zip

# These link to old harp csv files; they are now generated explicitly from raw data, below
# singularity exec ${UTILS_SIF} wget -O priors.tar.gz "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118196&authkey=AARukAbSz23Cz5w"
# singularity exec ${UTILS_SIF} tar -zxvf priors.tar.gz && rm priors.tar.gz 

# Download reference genome
if [ ! -f "${REFGENOME}" ]; then
    if [ ! -f "${REFGENOME_ZIPPED}" ]; then
        wget ${REFGENOME_URL}
    fi
    singularity exec ${UTILS_SIF} gunzip ${REFGENOME_ZIPPED}
else
    echo "Reference genome ${REFGENOME_FILENAME%.gz} already downloaded..."
fi

# If incorrect vcf md5 sum, remove corrupt file
if [ -f ${VCF_ZIPPED} ] && [ ! $(md5sum ${VCF_ZIPPED} | awk '{print $1}') == ${VCF_BGZIP_MD5} ]; then
    echo "removing improperly bgzipped VCF..."
    rm "${VCF_ZIPPED}" "${VCF}"
fi
if [ ! -f ${VCF_ZIPPED} ]; then
    echo "Downloading original VCF..."
    [ -f "${VCF}" ] && rm ${VCF}
    wget ${VCF_URL} && \
    if [ "${VCF_URL: -3}" == ".gz" ]; then
        gunzip ${VCF_ZIPPED}
        singularity exec ${UTILS_SIF} bgzip ${VCF}
    fi
else
    echo "Correct vcf exists and is bgzipped..."
fi

# Download repetitive region masking files

if [ ! -f repetitive.list ]; then
    mkdir -p tmp_repetitive/ && \
    wget ${REP_MASK_1}
    tar -zxvf chromOut.tar.gz -C tmp_repetitive/ && \
    rm chromOut.tar.gz

    wget ${REP_MASK_2}
    tar -zxvf chromTrf.tar.gz -C tmp_repetitive/ && \
    rm chromTrf.tar.gz
else
    echo "Repetitive masking file 'repetitive.list' already exists..."
fi

