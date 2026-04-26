#!/bin/bash
#SBATCH --job-name=taxa_clssification
#SBATCH --partition=compute
#SBATCH --ntasks-per-node=112
#SBATCH --nodes=1

# Create output directories
mkdir -p ./kraken_2out_class
mkdir -p ./kraken_2out_class/bracken

INPUT_DIR="./"
OUTPUT_DIR="./kraken_2out_class"
BRACKEN_OUTPUT_DIR="./kraken_2out_class/bracken"
DB_PATH="/hpcigib/raiyan.a/dengue_raw_fastq/k2_pluspfp/"


for prefix in $(ls ${INPUT_DIR}*_unmapped.1.fastq.gz | sed 's/_unmapped.1.fastq.gz//')
do
    # Run kraken2 
    kraken2 --db ${DB_PATH} \
    --threads 110 \
    --report ${OUTPUT_DIR}/$(basename ${prefix})_report.txt \
    --output ${OUTPUT_DIR}/$(basename ${prefix})_output.txt \
    --gzip-compressed --paired ${prefix}_unmapped.1.fastq.gz ${prefix}_unmapped.2.fastq.gz

    # Run bracken on kraken2 reports
    bracken -d ${DB_PATH} \
    -i ${OUTPUT_DIR}/$(basename ${prefix})_report.txt \
    -o ${BRACKEN_OUTPUT_DIR}/$(basename ${prefix})_bracken.txt
done

