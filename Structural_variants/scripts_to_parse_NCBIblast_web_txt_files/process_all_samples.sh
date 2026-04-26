#!/bin/bash
# Master script to process all BLAST alignment files
# CORRECT VERSION - handles multiple reads per sample properly
# Usage: ./process_all_samples.sh

set -e

echo "========================================================"
echo "BLAST Parser - Processing All Samples"
echo "========================================================"
echo ""

# Create output directory
mkdir -p blast_results
cd blast_results

# Copy scripts if needed
if [ ! -f "../parse_blast_correct.py" ]; then
    echo "Error: parse_blast_correct.py not found!"
    exit 1
fi

# Count files
total=$(ls ../*.txt 2>/dev/null | wc -l)
if [ $total -eq 0 ]; then
    echo "Error: No .txt files found in parent directory!"
    exit 1
fi

current=0
total_reads=0
total_samples=0

echo "Found $total BLAST alignment files to process"
echo ""

# Process each .txt file
for file in ../*.txt; do
    current=$((current + 1))
    basename=$(basename "$file" .txt)
    
    echo "[$current/$total] Processing: $basename"
    
    # Run the parser
    output=$(python3 ../parse_blast_correct.py "$file" "$basename" 2>&1)
    
    if [ $? -eq 0 ]; then
        # Extract number of reads from output
        reads=$(echo "$output" | grep "Total reads" | awk '{print $NF}')
        if [ -n "$reads" ]; then
            total_reads=$((total_reads + reads))
        fi
        total_samples=$((total_samples + 1))
        echo "  ✓ Parsed $reads reads"
    else
        echo "  ✗ Error processing $basename"
        echo "$output"
    fi
done

echo ""
echo "========================================================"
echo "Combining All Top1 Results with SampleName Column"
echo "========================================================"
echo ""

# Combine all top1 files
if [ -f "../combine_top1.py" ]; then
    python3 ../combine_top1.py combined_all_top1_hits.tsv
else
    echo "Warning: combine_top1.py not found, using awk method..."
    
    # Header with SampleName
    echo -e "SampleName\t$(head -1 $(ls *_top1.tsv | head -1))" > combined_all_top1_hits.tsv
    
    # Add data with sample names
    for file in *_top1.tsv; do
        sample=$(basename "$file" _top1.tsv)
        tail -n +2 "$file" | awk -v sample="$sample" '{print sample"\t"$0}' >> combined_all_top1_hits.tsv
    done
    
    combined_rows=$(tail -n +2 combined_all_top1_hits.tsv | wc -l)
    echo "✓ Combined $combined_rows rows"
fi

echo ""
echo "========================================================"
echo "Summary"
echo "========================================================"
echo ""
echo "Samples processed: $total_samples"
echo "Total reads across all samples: $total_reads"
echo ""

echo "Individual files per sample:"
ls -1 *_all.tsv 2>/dev/null | wc -l | xargs echo "  - All hits files:"
ls -1 *_top5.tsv 2>/dev/null | wc -l | xargs echo "  - Top5 files:"
ls -1 *_top1.tsv 2>/dev/null | wc -l | xargs echo "  - Top1 files:"

echo ""
echo "Combined file:"
if [ -f "combined_all_top1_hits.tsv" ]; then
    total_lines=$(wc -l < combined_all_top1_hits.tsv)
    data_lines=$((total_lines - 1))
    echo "  - combined_all_top1_hits.tsv: $data_lines rows"
    echo ""
    echo "Preview (first 5 rows):"
    head -5 combined_all_top1_hits.tsv | cut -f1-6
fi

echo ""
echo "========================================================"
echo "✓ All processing complete!"
echo "========================================================"
echo ""
echo "Results location: $(pwd)"
echo ""
echo "Next steps:"
echo "  1. Check: head combined_all_top1_hits.tsv"
echo "  2. Species count: tail -n +2 combined_all_top1_hits.tsv | cut -f4 | sort | uniq -c | sort -rn | head -20"
echo "  3. Reads per sample: tail -n +2 combined_all_top1_hits.tsv | cut -f1 | sort | uniq -c"
