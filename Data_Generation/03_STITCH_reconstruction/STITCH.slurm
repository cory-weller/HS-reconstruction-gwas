#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=4
#SBATCH --mem 64G
#SBATCH -t 0-16:00:00
#SBATCH -p largemem
#SBATCH --account berglandlab

module load gcc
module load R/3.5.1
module load htslib


population_sizes=(0 100 1000 5000)
N=${SLURM_ARRAY_TASK_ID}

population_name=${1}
population_size=${population_sizes[$N]}
n_generations=$(echo $population_name | cut -d "_" -f 4 | cut -d "F" -f 2)

# Create iteration folder
mkdir -p ${population_name}_${population_size} && cd ${population_name}_${population_size}

# Create bamlist
ls ../../02_simulate_reads/${population_name}/ | \
  grep ".bam$" | \
  sort -t "." -nk1,1 | \
  awk -F "." -v N=${population_size} '$1 <= N {print}' | \
  sed "s@^@../../02_simulate_reads/${population_name}/@g" > bamlist.txt

# Create sampleNames_file
cat bamlist.txt | sed 's@.*/@@g' > samples.txt

# Link to posfile
posfile="2L.sites"

# 2L.sites is a plain text file containing ref and alt alleles of dmel chromosome 2L, of the format
# 2L	4998	G	A
# 2L	5039	C	T
# 2L	5002	G	T
# 2L	5040	G	A


../STITCH/STITCH.R \
--inputBundleBlockSize=100 \
--method=pseudoHaploid \
--sampleNames_file=samples.txt \
--outputdir=./ \
--chr=2L \
--posfile=${posfile} \
--bamlist=bamlist.txt \
--K=32 \
--nCores=4 \
--nGen=${n_generations}
