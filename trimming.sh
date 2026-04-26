#!/bin/bash
#SBATCH --job-name=TRIM
#SBATCH --partition=compute
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1 

# Path to the adapter file
ADAPTERS="/home/raiyan.ali/.conda/envs/amr_rna/share/trimmomatic-0.39-2/adapters/custom_adap.fa"

OUTPUT_DIR="./trimming"
mkdir -p $OUTPUT_DIR

for R1 in *_R1_001.fastq.gz; do
    # Extract the base name (without _R1.fastq.gz)
    BASE=$(basename $R1 _R1_001.fastq.gz)

    # Define the corresponding R2 file
    R2="${BASE}_R2_001.fastq.gz"

    # Define output file names within the output directory
    R1_PAIRED="${OUTPUT_DIR}/${BASE}_R1_paired.fastq.gz"
    R1_UNPAIRED="${OUTPUT_DIR}/${BASE}_R1_unpaired.fastq.gz"
    R2_PAIRED="${OUTPUT_DIR}/${BASE}_R2_paired.fastq.gz"
    R2_UNPAIRED="${OUTPUT_DIR}/${BASE}_R2_unpaired.fastq.gz"

    # Log and summary file names within the output directory
    LOG="${OUTPUT_DIR}/${BASE}_trimmomatic.log"
    SUMMARY="${OUTPUT_DIR}/${BASE}_summary.txt"

    # Run Trimmomatic
    trimmomatic PE -threads $THREADS -summary $SUMMARY -validatePairs \
        $R1 $R2 \
        $R1_PAIRED $R1_UNPAIRED \
        $R2_PAIRED $R2_UNPAIRED \
        ILLUMINACLIP:$ADAPTERS:2:30:10  SLIDINGWINDOW:4:15 MINLEN:35
done

