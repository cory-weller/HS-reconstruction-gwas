#!/usr/bin/env bash

while getopts p:f:l:c: option; do
        case "${option}"
        in
        p) IFS=',' populations=($OPTARG);;
        f) first_ind=${OPTARG};;
        l) last_ind=${OPTARG};;
        c) chromosome=${OPTARG};;
        esac
done

echo for each population: ${populations[@]}
echo start with ${first_ind} and end with ${last_ind}
echo on chromosome ${chromosome}

module purge
module load gcc
module load R/3.5.1

for population in ${populations[@]}; do
  mkdir -p ${population}
  # create wide .haps file for each individual
  Rscript convert_haps_to_wide.R ../01_forward_simulator/${population}.haps ${population} ${first_ind} ${last_ind} ${chromosome}
done
