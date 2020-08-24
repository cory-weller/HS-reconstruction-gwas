#!/usr/bin/env bash

wget -O 73.test_run.bam "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118194&authkey=AJNA5CF5UX1bsLI"
wget -O 77.test_run.bam "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118191&authkey=AENQGxjwtHH3BD0"
wget -O 78.test_run.bam "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118190&authkey=AHNAMKodHtgt0e0"
wget -O 79.test_run.bam "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118192&authkey=ABxLTZCWwwuyjCY"
wget -O 81.test_run.bam "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118195&authkey=AFdak4Jwyei0NzE"
wget -O 82.test_run.bam "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118193&authkey=APj6cjxPbLdgl6c"

wget -O priors.tar.gz "https://onedrive.live.com/download?cid=77DD71E598E5B51B&resid=77DD71E598E5B51B%2118196&authkey=AARukAbSz23Cz5w"
tar -zxvf priors.tar.gz && rm priors.tar.gz 

# Prepare Reference

# Prepare PICARD sequence dictionary
module load picard/2.20.6

java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary \
  R=dgrp2.reference.fasta \
  O=dgrp2.reference.dict

# Prepare SAMTOOLS index
module load samtools/1.9
samtools faidx dgrp2.reference.fasta

# Prepare BWA index
module load bwa/0.7.17
bwa index dgrp2.reference.fasta

# tabix index vcf file
# if necessary, first bgzip vcf file
module load tabix/0.2.6
tabix -p vcf haplotypes.vcf.gz
tabix -p vcf haplotypes.polarized.vcf.gz


# split fasta file into separate files
bash splitFasta.sh dgrp2.reference.fasta
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

  buildindex(basename="dgrp2", reference="dgrp2.reference.fasta")

INNER
EOF