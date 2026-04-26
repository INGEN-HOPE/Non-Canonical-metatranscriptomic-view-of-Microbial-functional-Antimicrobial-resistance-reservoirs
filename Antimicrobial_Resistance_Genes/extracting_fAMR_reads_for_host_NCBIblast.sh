#!/bin/bash

set -euo pipefail

OUTDIR="extracted_fasta"

for GROUP in severe non-severe negative; do
  UNMAPPED_DIR="./${GROUP}/mapping/unmapped"

  echo "Processing group: $GROUP"

  for READS_FILE in *.reads.txt; do
    SAMPLE=$(basename "$READS_FILE" .reads.txt)

    FQ1="${UNMAPPED_DIR}/${SAMPLE}_unmapped.1.fastq.gz"
    FQ2="${UNMAPPED_DIR}/${SAMPLE}_unmapped.2.fastq.gz"

    # skip if fastq files do not exist in this group
    if [[ ! -f "$FQ1" || ! -f "$FQ2" ]]; then
      continue
    fi

    echo "  Extracting $SAMPLE from $GROUP"

    seqkit grep -n -f "$READS_FILE" "$FQ1" "$FQ2" \
      | seqkit fq2fa \
      -o "${OUTDIR}/${SAMPLE}.fasta"
  done
done

