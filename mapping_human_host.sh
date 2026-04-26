#!/bin/bash
#SBATCH --job-name=map
#SBATCH --partition=compute
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1 


# Path to HISAT2 index
HISAT2_INDEX="/lustre/raiyan.ali/hisat2_index/hisat2_index"

mkdir -p ./mapping
mkdir -p ./mapping/unmapped

# Loop through all R1 paired FASTQ files
for r1 in *_R1_paired.fastq.gz; do
    # Determine the corresponding R2 file
    r2="${r1/_R1_paired.fastq.gz/_R2_paired.fastq.gz}"

    # Check if the R2 file exists
    if [[ -f "$r2" ]]; then
        # Extract the base filename
        base=$(basename "$r1" _R1_paired.fastq.gz)

        # Run HISAT2
        hisat2 -p $THREADS -x $HISAT2_INDEX -1 "$r1" -2 "$r2" \
            -S "./mapping/${base}.sam" \
            --un-conc "./mapping/unmapped/${base}_unmapped.fastq" \
            --summary-file "./mapping/${base}_align_summary.txt"

        # Convert SAM to BAM, and sort BAM file
        samtools view -bS "./mapping/${base}.sam" | samtools sort -o "./mapping/${base}.sorted.bam"

        # Gzip the unpaired FASTQ files
        gzip "./mapping/unmapped/${base}_unmapped.1.fastq" "./mapping/unmapped/${base}_unmapped.2.fastq"

        # Clean up SAM file
        rm "./mapping/${base}.sam"
    else
        echo "Warning: R2 file for $r1 does not exist."
    fi
done

