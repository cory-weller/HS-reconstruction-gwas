#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem 8G
#SBATCH -t 0-12:00:00
#SBATCH -p standard
#SBATCH --account berglandlab

population_stem=${1}
chromosome=${2}

module load tabix

generation=${SLURM_ARRAY_TASK_ID}

population=${population_stem}_F${generation}

haplotypes_gzfile="../input_data/haplotypes.vcf.gz"

# Set up associative array of lineID : column, to be used by TABIX
vals=($(zgrep -m 1 "^#CHROM" ${haplotypes_gzfile} | xargs | cut -f 10-))
declare -A founderIndices

for i in $(seq 0 "${#vals[@]}"); do
    let j=$i+1
    founderIndices[[${vals[$i]}]]=$j
done

cd ${population}

for filename in *.haps.wide; do
  individual_n=${filename%.haps.wide}
  if [[ -f  ${individual_n}.geno ]]; then
    echo ${individual_n} already done
  else
    echo ${individual_n}
    while read chromosomeN start stop par1 par2; do
        if [[ $chromosomeN == $chromosome ]]; then
        col1=${founderIndices[[$par1]]}
        col2=${founderIndices[[$par2]]}
        tabix ../${haplotypes_gzfile} ${chromosomeN}:${start}-${stop} | awk -v col1=$col1 -v col2=$col2 '{print $1,$2,$col1,$col2}' >> ${individual_n}.geno
        fi
    done < <(awk 'NR > 1 {print}' ${filename})
  fi
done



# remove inside 'if' block to instead do all chromosomes
