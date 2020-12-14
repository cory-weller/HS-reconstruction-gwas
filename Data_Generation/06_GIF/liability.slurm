#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH --ntasks-per-node=5
#SBATCH --mem 50G
#SBATCH -t 0-10:00:00
#SBATCH -p standard
#SBATCH --account berglandlab

module load gcc
module load openmpi
module load R/3.5.1
module load tabix

slurm_array_id=${SLURM_ARRAY_TASK_ID}

population=${1}
replicates=${2}
n_loci=${3}
effect_size=${4}
maf_min_filter=${5}

#folder_name="${population}_${count}_${n_loci}_${effect_size}"

#mkdir -p data/ && cd data/

echo "running liability.GIF.R ${population} ${replicates} ${n_loci} ${effect_size} ${maf_min_filter}"

Rscript liability.GIF.R ${population} ${replicates} ${n_loci} ${effect_size} ${slurm_array_id} ${maf_min_filter}

# many loci small effect
# sbatch --array=1-10%1 liability.slurm inbredLines_128_F0 5 10 0.2 # DGRP 8244500
# sbatch --array=1-10%1 liability.slurm inbredLines_32_F0 5 10 0.2 # DGRPsmall 8271297
# sbatch --array=1-10%1 liability.slurm hybrid_swarm_32_F5 5 10 0.2 # HS_32 8244501
# sbatch --array=1-10%1 liability.slurm hybrid_swarm_128_F5 5 10 0.2 #  HS_128 8244502
# sbatch --array=1-10%1 liability.slurm hybrid_swarm_128_F1 5 10 0.2  #  HS_128 8271304
# sbatch --array=1-10%1 liability.slurm hybrid_swarm_128_F2 5 10 0.2 #  HS_128  8271305
# sbatch --array=1-10%1 liability.slurm RILs_8_F50 5 10 0.2  # DSPR 8244503
# sbatch --array=1-10%2 liability.slurm outbred_128_F50 5 10 0.2 # Outbred 8248166
# sbatch --array=1-10%1 liability.slurm hybrid_swarm_32_F2 5 10 0.2 # HS_F2 8244505
# sbatch --array=1-10%1 liability.slurm hybrid_swarm_32_F1 5 10 0.2  # HS_F1 8244506

# Single locus large effect
# sbatch --array=1-1%1 liability.slurm hybrid_swarm_128_F0 5 1 0.4 0.05 # DGRP 0.05 MAF
# sbatch --array=1-1%1 liability.slurm hybrid_swarm_128_F0 5 1 0.4 0 # DGRP all sites



# sbatch --array=1-10%1 liability.slurm inbredLines_128_F0 5 1 0.4 # DGRP 8343358
# sbatch --array=1-10%1 liability.slurm inbredLines_32_F0 5 1 0.4 # DGRPsmall 8343359
# sbatch --array=1-10%1 liability.slurm hybrid_swarm_32_F5 5 1 0.4 # HS_32 8343360
# sbatch --array=1-10%1 liability.slurm hybrid_swarm_128_F5 5 1 0.4 #  HS_128 8343361
# sbatch --array=1-10%1 liability.slurm hybrid_swarm_128_F1 5 1 0.4  #  HS_128 8343362
# sbatch --array=1-10%1 liability.slurm hybrid_swarm_128_F2 5 1 0.4 #  HS_128  8343363
# sbatch --array=1-10%1 liability.slurm RILs_8_F50 5 1 0.4  # DSPR 8343364
# sbatch --array=1-10%1 liability.slurm outbred_128_F50 5 1 0.4 # Outbred 8343365
# sbatch --array=1-10%1 liability.slurm hybrid_swarm_32_F2 5 1 0.4 # HS_F2 8343366
# sbatch --array=1-10%1 liability.slurm hybrid_swarm_32_F1 5 1 0.4  # HS_F1 8343367

# # rerun bad groups?
# sbatch --array=4-10%1 liability.slurm hybrid_swarm_128_F1 5 1 0.4  # 8375598
# sbatch --array=4 liability.slurm hybrid_swarm_32_F1 5 1 0.4  # 8375599
# sbatch --array=5-6%1 liability.slurm hybrid_swarm_128_F2 5 1 0.4 #  8375600
# sbatch --array=5-6%1 liability.slurm inbredLines_128_F0 5 1 0.4 #  8375601
# sbatch --array=4-5%1 liability.slurm inbredLines_32_F0 5 1 0.4 # 8375604
# sbatch --array=6-10%1 liability.slurm outbred_128_F50 5 1 0.4
