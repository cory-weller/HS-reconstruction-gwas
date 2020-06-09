#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 8G
#SBATCH -t 0-12:00:00
#SBATCH -p standard
#SBATCH --account berglandlab

# This script converts haplotype map format to site-by-site genotypes

module load tabix

population_stem="${1}"
chromosome="${2}"
first_ind=${3}
last_ind=${4}
generation="${SLURM_ARRAY_TASK_ID}"

population="${population_stem}_F${generation}"

haplotypes_gzfile="../input_data/haplotypes.vcf.gz"

mkdir -p "${population}/imputed_genotypes/"

for ind_id in $(seq ${first_ind} ${last_ind}); do

haplotypes_in_file="${population}/imputed_haplotypes/${ind_id}.${chromosome}.estimate.haps"
genotypes_out_file="${population}/imputed_genotypes/${ind_id}.${chromosome}.estimate.geno"


# Set up associative array of lineID : column, to be used by TABIX
vals=($(zgrep -m 1 "^#CHROM" ${haplotypes_gzfile} | xargs | cut -f 10-))
declare -A founderIndices

for i in $(seq 0 "${#vals[@]}"); do
    let j=$i+1
    founderIndices[[${vals[$i]}]]=$j
done

> ${genotypes_out_file}

while read chromosomeN start stop par1 par2; do
    if [[ $chromosomeN == $chromosome ]]; then
    col1=${founderIndices[[$par1]]}
    col2=${founderIndices[[$par2]]}
    tabix ${haplotypes_gzfile} ${chromosomeN}:${start}-${stop} | awk -v col1=$col1 -v col2=$col2 '{print $1,$2,$col1,$col2}' >> ${genotypes_out_file}
    fi
done < <(awk 'NR > 1 {print}' ${haplotypes_in_file})

done
