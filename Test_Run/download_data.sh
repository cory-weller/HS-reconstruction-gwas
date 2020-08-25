#!/usr/bin/env bash

UTILS_SIF="../utils.sif"
HARP_SIF="../harp.sif"
VCF_URL="http://dgrp2.gnets.ncsu.edu/data/website/dgrp2.vcf"
VCF_FILENAME=$(singularity exec ${UTILS_SIF} basename ${VCF_URL})
REFGENOME_URL="ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/dmel_r5.57_FB2014_03/fasta/dmel-all-chromosome-r5.57.fasta.gz"

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
REFGENOME_FILENAME=$(singularity exec ${UTILS_SIF} basename ${REFGENOME_URL})
if [[ ${REFGENOME_FILENAME: -3} == ".gz" ]]; then
  singularity exec ${UTILS_SIF} gunzip ${REFGENOME_FILENAME}
  REFGENOME_FILENAME=${REFGENOME_FILENAME%.gz}
fi

# Create reference sequence dictionary
singularity exec ${UTILS_SIF} java -jar /opt/gatk.jar CreateSequenceDictionary \
  -R "${REFGENOME_FILENAME}" \
  -O "${REFGENOME_FILENAME%.fasta}".dict

# Prepare SAMTOOLS index
singularity exec ${HARP_SIF} samtools faidx ${REFGENOME_FILENAME}

# Prepare BWA index
singularity exec ${UTILS_SIF} bwa index ${REFGENOME_FILENAME}

# Download DGRP freeze 2 VCF
singularity exec ${UTILS_SIF} wget ${VCF_URL}
singularity exec ${UTILS_SIF} bgzip ${VCF_FILENAME}
BGZIP_VCF="${VCF_FILENAME}.gz"

# Sort VCF
singularity exec ${UTILS_SIF} java -Xmx2G -jar /opt/gatk4.jar SortVcf \
    -I ${BGZIP_VCF} \
    -O ${BGZIP_VCF%.vcf.gz}.sorted.vcf \
    --SEQUENCE_DICTIONARY ${REFGENOME_FILENAME%.fasta}.dict \
    --MAX_RECORDS_IN_RAM 100000

# Prepare VCF subset without indels
singularity exec ${UTILS_SIF} awk -v file=${BGZIP_VCF%.vcf.gz}.sorted.noIndel.vcf -F "\t" '{
if($0 ~ /^#/ || $3 ~ /_SNP/)
    print >> file}' < ${BGZIP_VCF%.vcf.gz}.sorted.vcf

# Prepare VCF subset with only indels
singularity exec ${UTILS_SIF} awk -v file=${BGZIP_VCF%.vcf.gz}.sorted.indel.vcf -F "\t" '{
if($0 ~ /^#/ || $3 ~ /_DEL/ || $3 ~ /_INS/)
    print >> file}' < ${BGZIP_VCF%.vcf.gz}.sorted.vcf

### STOP HERE FOR NOW




























# bgzip subset VCF
singularity exec ${UTILS_SIF} bgzip ${BGZIP_VCF%.vcf.gz}.noIndel.vcf
singularity exec ${UTILS_SIF} bgzip ${BGZIP_VCF%.vcf.gz}.indel.vcf

# tabix index VCF
singularity exec ${UTILS_SIF} tabix -p vcf ${BGZIP_VCF}


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