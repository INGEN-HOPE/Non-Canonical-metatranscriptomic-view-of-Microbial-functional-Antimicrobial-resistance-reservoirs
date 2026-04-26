#!/bin/bash
#SBATCH --job-name=AR_G
#SBATCH --partition=smp
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1

for file in /lustre/raiyan.ali/sev_vs_non-sev/severe/mapping/unmapped/*_unmapped.1.fastq.gz; do
    # Strip directory and suffix
    base=$(basename "$file" _unmapped.1.fastq.gz)

    rgi bwt \
        -1 "/lustre/raiyan.ali/sev_vs_non-sev/severe/mapping/unmapped/${base}_unmapped.1.fastq.gz" \
        -2 "/lustre/raiyan.ali/sev_vs_non-sev/severe/mapping/unmapped/${base}_unmapped.2.fastq.gz" \
        --output_file "/lustre/raiyan.ali/sev_vs_non-sev/severe/card_output/${base}" \
        --local -n 40 --include_wildcard -a kma
done

