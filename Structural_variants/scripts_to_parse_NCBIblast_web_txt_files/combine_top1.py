#!/usr/bin/env python3
"""
Combine all *_top1.tsv files into one file with SampleName column
Updated to handle NA entries for reads with no hits
"""

import sys
import os
import glob

def combine_top1_files(output_file='combined_all_top1_hits.tsv'):
    """Combine all *_top1.tsv files with SampleName column added."""

    # Find all top1 files
    top1_files = sorted(glob.glob('*_top1.tsv'))

    if not top1_files:
        print("Error: No *_top1.tsv files found in current directory!")
        sys.exit(1)

    print(f"Found {len(top1_files)} top1 files to combine")

    # Read header from first file
    with open(top1_files[0], 'r') as f:
        header = f.readline().strip()

    # Write combined file with SampleName column
    with open(output_file, 'w') as outfile:
        # Write header with SampleName as first column
        outfile.write('SampleName\t' + header + '\n')

        total_rows = 0
        total_with_hits = 0
        total_no_hits = 0
        sample_stats = {}

        for top1_file in top1_files:
            # Extract sample name from filename (remove _top1.tsv)
            sample_name = top1_file.replace('_top1.tsv', '')
            
            sample_total = 0
            sample_hits = 0
            sample_no_hits = 0

            with open(top1_file, 'r') as infile:
                # Skip header
                infile.readline()

                # Add each row with sample name prepended
                for line in infile:
                    line = line.strip()
                    if line:
                        outfile.write(sample_name + '\t' + line + '\n')
                        total_rows += 1
                        sample_total += 1
                        
                        # Check if this is a hit or NA
                        fields = line.split('\t')
                        if len(fields) >= 3:
                            # Check ScientificName column (index 2 in original, index 3 with SampleName)
                            if fields[2] == 'NA':
                                total_no_hits += 1
                                sample_no_hits += 1
                            else:
                                total_with_hits += 1
                                sample_hits += 1
            
            # Store sample statistics
            sample_stats[sample_name] = {
                'total': sample_total,
                'hits': sample_hits,
                'no_hits': sample_no_hits
            }

        print(f"\n✓ Combined {total_rows} rows from {len(top1_files)} samples")
        print(f"✓ Output: {output_file}")
        print(f"\nOverall Statistics:")
        print(f"  - Total reads: {total_rows}")
        print(f"  - Reads with hits: {total_with_hits} ({100*total_with_hits/total_rows:.1f}%)" if total_rows > 0 else "")
        print(f"  - Reads with NO hits (NA): {total_no_hits} ({100*total_no_hits/total_rows:.1f}%)" if total_rows > 0 else "")
        
        # Show samples with high no-hit rates (>50%)
        high_no_hit_samples = []
        for sample, stats in sample_stats.items():
            if stats['total'] > 0:
                no_hit_pct = 100 * stats['no_hits'] / stats['total']
                if no_hit_pct > 50:
                    high_no_hit_samples.append((sample, no_hit_pct, stats))
        
        if high_no_hit_samples:
            print(f"\n⚠️  Samples with >50% no-hit reads:")
            for sample, pct, stats in sorted(high_no_hit_samples, key=lambda x: x[1], reverse=True)[:10]:
                print(f"  - {sample}: {pct:.1f}% ({stats['no_hits']}/{stats['total']})")
        
        print(f"\nColumns: SampleName, ReadName, Description, ScientificName, CommonName, ...")

if __name__ == '__main__':
    output_file = sys.argv[1] if len(sys.argv) > 1 else 'combined_all_top1_hits.tsv'
    combine_top1_files(output_file)
