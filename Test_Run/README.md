#

## Preparing data
1. Ensure the containers `harp.simg` and `utils.simg` are downloaded in the previous directory.

```bash
./download_data.sh  # retrieves a handful of necessary files
./prep_data.sh  # builds reference sequence dict, samtools index, sorts VCF file, etc.

sample_ids=( *_R1_001.fastq )
sample_ids=( ${filenames[@]%_R1_001.fastq} )

for sample in sample_ids; do
    ./map_reads.sh ${sample}
done

# batch reconstruct all FILENAME.final.bam in job array with SLURM
./batch_reconstruct.sh

# or submit one at a time to SLURM
source PARAMETERS.config && sbatch ${SLURM_OPTS} ./reconstruct.sh FILENAME.final.bam

# or run one job locally
./reconstruct.sh FILENAME.final.bam
```