#!/usr/bin/env bash

UTILS_SIF="../utils.sif"
HARP_SIF="../harp.sif"


# VCF_URL="http://dgrp2.gnets.ncsu.edu/data/website/dgrp2.vcf" # gzipped version from zenodo faster to download (equivalent md5 sum)
VCF_URL="https://zenodo.org/api/files/eef9cec5-1bf4-498e-8ac6-e90aa8c8ab1c/dgrp2.vcf.gz"
VCF_MD5="af831816f41fe4a779e9abb44c5cdc62"
VCF_BGZIP_MD5="f4d8551dacf10b6e790d53deeed1f02a"
VCF_GZIP_MD5="a9ece6b8c4c6b8aaf6797874eaddb369"
VCF_FILENAME=$(singularity exec ${UTILS_SIF} basename ${VCF_URL})
VCF_SORTED_MD5=""

# Reference Genome: Drosophila melanogaster release 5, version 57
REFGENOME_URL="ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/fasta/dmel-all-chromosome-r5.57.fasta.gz"
REFGENOME_FILENAME=$(singularity exec ${UTILS_SIF} basename ${REFGENOME_URL})
REFGENOME_GZIP_MD5="f6e168a1fe8741be2fdce6c2d0602b41"
REFGENOME_MD5="203ef91d30c76497cd1391a0d2f92827"

# download bam files
# wget -O 73.test_run.bam "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118194&authkey=AJNA5CF5UX1bsLI"
# wget -O 77.test_run.bam "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118191&authkey=AENQGxjwtHH3BD0"
# wget -O 78.test_run.bam "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118190&authkey=AHNAMKodHtgt0e0"
# wget -O 79.test_run.bam "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118192&authkey=ABxLTZCWwwuyjCY"
# wget -O 81.test_run.bam "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118195&authkey=AFdak4Jwyei0NzE"
# wget -O 82.test_run.bam "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118193&authkey=APj6cjxPbLdgl6c"

# Download low coverage reads
singularity exec ${UTILS_SIF} wget -O low_coverage_reads.zip "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118197&authkey=AAJg8Dw7UluE1w8" && \
singularity exec ${UTILS_SIF} unzip low_coverage_reads.zip && \
rm low_coverage_reads.zip

# TODO change this to conversion from dgrp2 vcf file
singularity exec ${UTILS_SIF} wget -O priors.tar.gz "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118196&authkey=AARukAbSz23Cz5w"
singularity exec ${UTILS_SIF} tar -zxvf priors.tar.gz && rm priors.tar.gz 

# Download reference genome
singularity exec ${UTILS_SIF} wget ${REFGENOME_URL}
if [[ ${REFGENOME_FILENAME: -3} == ".gz" ]]; then
  singularity exec ${UTILS_SIF} gunzip ${REFGENOME_FILENAME}
  REFGENOME_FILENAME=${REFGENOME_FILENAME%.gz}
fi

# Create reference sequence dictionary
singularity exec ${UTILS_SIF} java -jar /opt/gatk4.jar CreateSequenceDictionary \
  -R "${REFGENOME_FILENAME}" \
  -O "${REFGENOME_FILENAME%.fasta}".dict

# Prepare SAMTOOLS index
singularity exec ${HARP_SIF} samtools faidx ${REFGENOME_FILENAME}

# Prepare BWA index
singularity exec ${UTILS_SIF} bwa index ${REFGENOME_FILENAME}

# Download DGRP freeze 2 VCF
while true; do
  if [[ -f ${VCF_FILENAME} ]]; then
    if [[ $(md5sum ${VCF_FILENAME} | awk '{print $1}') == ${VCF_BGZIP_MD5} ]]; then
      echo "BGZIP VCF is ready"
      break
    elif [[ ! $(md5sum ${VCF_FILENAME} | awk '{print $1}') == ${VCF_GZIP_MD5} ]]; then
      echo "Removing ${VCF_FILENAME} due to corrupted md5 checksum"
      rm ${VCF_FILENAME}
      continue
    elif [[ $(md5sum ${VCF_FILENAME} | awk '{print $1}') == ${VCF_GZIP_MD5} ]]; then
      echo "Unzipping non-BGZIP file..."
      gunzip ${VCF_FILENAME}
    fi
  elif [[ -f ${VCF_FILENAME%.gz} ]]; then
    if [[ $(md5sum ${VCF_FILENAME%.gz} | awk '{print $1}') == ${VCF_MD5} ]]; then
      echo "Compressing with BGZIP..."
      singularity exec ${UTILS_SIF} bgzip ${VCF_FILENAME%.gz}
    else
      echo "Removing ${VCF_FILENAME%.gz} due to corrupted md5 checksum"
      rm ${VCF_FILENAME%.gz}
      continue
    fi
    continue
  else
    echo "Downloading ${VCF_FILENAME}"
    singularity exec ${UTILS_SIF} wget ${VCF_URL}
  fi
done

# Sort VCF

if [[ ! -f ${VCF_FILENAME%.vcf.gz}.sorted.vcf ]]; then
  echo "Sorting VCF..."
  # dgrp2.sorted.vcf md5sum: dac0ddb8da32d9e4b68458b7aff28caf
  singularity exec ${UTILS_SIF} java -Xmx2G -jar /opt/gatk4.jar SortVcf \
      -I ${VCF_FILENAME} \
      -O ${VCF_FILENAME%.vcf.gz}.sorted.vcf \
      --SEQUENCE_DICTIONARY ${REFGENOME_FILENAME%.fasta}.dict \
      --MAX_RECORDS_IN_RAM 100000
  echo "BGZIP on sorted VCF..."
  singularity exec ${UTILS_SIF} bgzip ${VCF_FILENAME%.vcf.gz}.sorted.vcf
  # dgrp2.sorted.vcf.gz md5sum: e8364935cb63d5fc945a072321d8631c
fi

# Prepare VCF subset without indels
singularity exec ${UTILS_SIF} awk -v file=${VCF_FILENAME%.vcf.gz}.sorted.noIndel.vcf -F "\t" '{
if($0 ~ /^#/ || $3 ~ /_SNP/)
    print >> file}' < <(zcat ${VCF_FILENAME%.vcf.gz}.sorted.vcf.gz)


# Prepare VCF subset with only indels
singularity exec ${UTILS_SIF} awk -v file=${VCF_FILENAME%.vcf.gz}.sorted.indel.vcf -F "\t" '{
if($0 ~ /^#/ || $3 ~ /_DEL/ || $3 ~ /_INS/)
    print >> file}' < <(zcat ${VCF_FILENAME%.vcf.gz}.sorted.vcf.gz)

exit 0
### STOP HERE FOR NOW
# dgrp2.sorted.indel.vcf md5: 0fe9296e7c1681697492c9175a43986f
# dgrp2.sorted.noIndel.vcf md5: c4793c71cf076964441944845ab8b921


# Download repetitive region masking files
mkdir -p tmp_repetitive/ && \
singularity exec ${UTILS_SIF} wget http://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/chromOut.tar.gz
tar -zxvf chromOut.tar.gz -C tmp_repetitive/ && \
rm chromOut.tar.gz

singularity exec ${UTILS_SIF} wget http://hgdownload.soe.ucsc.edu/goldenPath/dm3/bigZips/chromTrf.tar.gz
tar -zxvf chromTrf.tar.gz -C tmp_repetitive/ && \
rm chromTrf.tar.gz

# Format as samtools ranges, e.g. CHROM:START-STOP
# repeatmasker is 1-based, closed [start, end] 
# bed is  0-based, half-open [start-1, end)
# Format as 1-based .intervals and concatenate
for CHR in 2L 2R 3L 3R X; do
  cat <(awk 'BEGIN{OFS="\t"} NR > 3 {printf ("%s:%s-%s\n",$5,$6,$7)}' tmp_repetitive/${CHR}/chr${CHR}.fa.out) \
      <(awk 'BEGIN{OFS="\t"} {printf ("%s:%s-%s\n", $1,$2+1,$3)}' tmp_repetitive/trfMaskChrom/chr${CHR}.bed) | \
      sort -nk 2,3 | \
      sed "s/chr//g" >> repetitive.list
done && \
rm -rf tmp_repetitive/

# Index sorted vcf
singularity exec ${UTILS_SIF} java -jar /opt/gatk4.jar IndexFeatureFile \
  --input ${VCF_FILENAME%.vcf.gz}.sorted.noIndel.vcf

# Exclude repeptitive intervals
singularity exec ${UTILS_SIF} java -Xmx2G -jar /opt/gatk4.jar SelectVariants \
  --variant  ${VCF_FILENAME%.vcf.gz}.sorted.noIndel.vcf \
  --exclude-intervals repetitive.list \
  --output ${VCF_FILENAME%.vcf.gz}.sorted.noIndel.noRep.vcf


# Generate HARP genotype csv
singularity exec ${UTILS_SIF} python3 generate_harp_csv.py ${VCF_FILENAME%.vcf.gz}.sorted.noIndel.noRep.vcf


# Generate associative array for chromosome lengths
awk '{print $1,$2}' dmel-all-chromosome-r5.57.fasta.fai


# Read chromosome lengths into associative array from the fasta index
declare -A LENGTHS
while read -r CHROM LENGTH; do
  LENGTHS[${CHROM}]=${LENGTH}
done < <(awk '{print $1,$2}' ${REFGENOME_FILENAME}.fai)


# harp parameters
HARP_STEP=100000
HARP_WIDTH=100000
chromosome="2L"
chrLength="23011544"
harp="/scratch/caw5cv/genome-reconstruction-revision/04_RABBIT/harp"
REFGENOME_FILENAME="${projectDir}/input_data/dgrp2.reference.fasta"
bamFolder="${projectDir}/02_simulate_reads/${population}/"
priorsGzFile="${projectDir}/04_RABBIT/${population}/${chromosome}.subset.priors.csv.gz"

BAM_FILENAME
CHROMOSOME
REFGENOME_FILENAME
HARP_CSV_FILENAME
ORIG_DIR=$(pwd)
function getHarpFreqs {
  BAM_FILENAME=${1}
  CHROMOSOME=${2}
  CHROMOSOME_LENGTH=${LENGTHS[${CHROMOSOME}]}

  CHROMOSOME=${3}
  tmpWorkDir=${4}
  outputDir=${5}
  chrLength=${6}
  REFGENOME_FILENAME=${7}
  priors=${8}
  window_step=${9}
  window_width=${10}
  harp=${11}

  HARP_CSV_FILENAME="${CHROMOSOME}.harp.csv"

  BAM="${ORIG_DIR}/${BAM_FILENAME}"
  REFGENOME="${ORIG_DIR}/${REFGENOME_FILENAME}"
  HARP_CSV="${ORIG_DIR}/${HARP_CSV_FILENAME}"

  # make tmp work directory
  echo creating directory "${tmpWorkDir}/${ind}"
  mkdir -p "${tmpWorkDir}/${ind}" && cd "${tmpWorkDir}/${ind}"

  # Run harp like
  echo running harp like
  singularity exec ${UTILS_SIF} /opt/harp like \
  --bam ${BAM} \
  --region ${CHROMOSOME}:1-${CHROMOSOME_LENGTH} \
  --refseq ${REFGENOME} \
  --snps ${HARP_CSV} \
  --stem $ind.$CHROMOSOME

  # Run harp freq
  echo running harp freq
 singularity exec ${UTILS_SIF} /opt/harp freq \
  --bam ${BAM} \
  --region ${CHROMOSOME}:1-${CHROMOSOME_LENGTH} \
  --refseq ${REFGENOME} \
  --snps ${HARP_CSV} \
  --stem $ind.$CHROMOSOME \
  --window_step ${HARP_STEP} \
  --window_width ${HARP_WIDTH} \
  --em_min_freq_cutoff 0.0001

  # Cleanup and move out
  echo cleaning up
  cp "${ind}.${CHROMOSOME}.freqs" ${outputDir} && \
  cd ${tmpWorkDir} && \
  rm -rf ${tmpWorkDir}/${ind}
  echo done
}

export -f getHarpFreqs


















# bgzip subset VCF
singularity exec ${UTILS_SIF} bgzip ${VCF_FILENAME%.vcf.gz}.noIndel.vcf
singularity exec ${UTILS_SIF} bgzip ${VCF_FILENAME%.vcf.gz}.indel.vcf

# tabix index VCF
singularity exec ${UTILS_SIF} tabix -p vcf ${VCF_FILENAME}


# split fasta file into separate files
bash splitFasta.sh ${REFGENOME_FILENAME}
# creates 2L.fa, 2R.fa, ... etc

# prepare {chromosome}.sites files
zcat haplotypes.vcf.gz | awk -F "\t" 'BEGIN{OFS="\t";} $1 !~ /^#/ {s=$1".sites"; print $1,$2,$4,$5 > s}'

# prepare heterozygous VCF file for 2L (for ASEReadCounter)
awk -F "\t"  '
BEGIN{OFS="\t";}
{
if ($1 ~ /^##/)
	print $0;
else if ($1 ~ /^#CHROM/)
	print $1,$2,$3,$4,$5,$6,$7,$8,$9,"HET";
else if ($1 ~ /^2L/)
	print $1,$2,$3,$4,$5,$6,$7,$8,$9,"0/1";
}' <(zcat haplotypes.vcf.gz) > variants_het_2L.vcf 

module load gatk
gatk IndexFeatureFile --feature-file variants_het_2L.vcf

# Create Rsubread Index
sbatch <<EOF
#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 30G
#SBATCH -t 0-0:15:00
#SBATCH -p largemem
#SBATCH --account berglandlab

Rscript - <<INNER

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("Rsubread")) {
  BiocManager::install("Rsubread")
  library(Rsubread)
}

  buildindex(basename="dgrp2", reference="${REFGENOME_FILENAME}")

INNER
EOF