#!/usr/bin/env python3
"""
Parse BLAST alignment output to TSV format with TOP HITS
CORRECTLY FIXED VERSION - NOW HANDLES "No significant similarity found"
- Each file = 1 sample with multiple reads (queries)
- Each read has 0-100 BLAST hits
- Queries with no hits get NA entries
- Generates: all hits, top 5 per read, top 1 per read
"""

import sys
import re
from collections import defaultdict

def parse_evalue(evalue_str):
    """Convert E-value string to float for comparison."""
    try:
        return float(evalue_str)
    except ValueError:
        return 1.0

def parse_query_cover(qcov_str):
    """Convert query coverage string (e.g., '84%') to float."""
    try:
        return float(qcov_str.rstrip('%'))
    except ValueError:
        return 0.0

def parse_perc_ident(pid_str):
    """Convert percent identity string to float."""
    try:
        return float(pid_str)
    except ValueError:
        return 0.0

def score_hit(hit):
    """
    Calculate composite score for ranking hits.
    Priority: PercIdent (highest), QueryCover (highest), Evalue (lowest)
    """
    perc_ident = parse_perc_ident(hit['PercIdent'])
    query_cover = parse_query_cover(hit['QueryCover'])
    evalue = parse_evalue(hit['Evalue'])

    # E-value: lower is better, convert to score (higher is better)
    if evalue == 0:
        evalue_score = 1000
    else:
        try:
            evalue_score = -1 * abs(float(f"{evalue:.2e}".split('e')[1]))
        except:
            evalue_score = 0

    # Weighted score: PercIdent most important, then QueryCover, then Evalue
    score = (perc_ident * 100) + (query_cover * 10) + (evalue_score * 1)
    return score

def is_valid_hit_line(fields):
    """Validate that a line is a proper BLAST hit line."""
    if len(fields) < 10:
        return False

    try:
        # Taxid should be numeric
        if not fields[-8].replace('-', '').isdigit():
            return False

        # Max score should be numeric
        if not fields[-7].replace('.', '').replace('-', '').isdigit():
            return False

        # Query cover should end with %
        if not fields[-5].endswith('%'):
            return False

        # Percent identity should be numeric
        if not fields[-3].replace('.', '').isdigit():
            return False

        return True
    except (IndexError, ValueError):
        return False

def create_na_hit(read_name):
    """Create a hit dictionary with NA values for queries with no similarity."""
    return {
        'ReadName': read_name,
        'Description': 'NA',
        'ScientificName': 'NA',
        'CommonName': 'NA',
        'Taxid': 'NA',
        'MaxScore': 'NA',
        'TotalScore': 'NA',
        'QueryCover': 'NA',
        'Evalue': 'NA',
        'PercIdent': 'NA',
        'AccLen': 'NA',
        'Accession': 'NA'
    }

def parse_blast_file(input_file, output_base):
    """
    Parse BLAST alignment file and write to three TSV files.

    For each sample file:
    - Reads all queries (reads) and their hits
    - Generates all hits file
    - Generates top 5 hits per read
    - Generates top 1 hit per read
    - Handles queries with no significant similarity (adds NA entries)
    """

    all_hits_file = output_base + '_all.tsv'
    top5_file = output_base + '_top5.tsv'
    top1_file = output_base + '_top1.tsv'

    # Store hits grouped by read name
    hits_by_read = defaultdict(list)
    total_queries = 0
    total_hits = 0
    no_hit_queries = 0

    with open(input_file, 'r') as infile:
        current_read = None
        skip_next = 0
        no_similarity_flag = False

        for line in infile:
            line = line.rstrip('\n')

            if not line.strip():
                continue

            # Detect new query
            if line.startswith('Query #'):
                # If previous query had no hits, ensure it's recorded
                if current_read and current_read not in hits_by_read:
                    hits_by_read[current_read] = []
                    no_hit_queries += 1
                
                match = re.search(r'Query #\d+: (.+?) Query ID:', line)
                if match:
                    current_read = match.group(1).strip()
                    total_queries += 1
                    skip_next = 2  # Skip next 2 lines (blank + header)
                    no_similarity_flag = False
                continue

            # Detect "No significant similarity found"
            if current_read and 'No significant similarity found' in line:
                no_similarity_flag = True
                if current_read not in hits_by_read:
                    hits_by_read[current_read] = []
                    no_hit_queries += 1
                continue

            # Skip header lines after Query
            if skip_next > 0:
                skip_next -= 1
                continue

            # Skip the Description/Scientific/Common header line
            if 'Scientific' in line and 'Common' in line and 'Taxid' in line:
                continue

            # Stop processing at alignment details
            if line.strip().startswith('>') or 'Alignments' in line:
                continue

            # Process data lines
            if current_read and len(line) > 60:
                fields = line.split()

                # Only process valid hit lines
                if not is_valid_hit_line(fields):
                    continue

                if len(fields) >= 10:
                    accession = fields[-1]
                    acc_len = fields[-2]
                    perc_ident = fields[-3]
                    evalue = fields[-4]
                    query_cover = fields[-5]
                    total_score = fields[-6]
                    max_score = fields[-7]
                    taxid = fields[-8]
                    common_name = fields[-9]

                    # Determine description end
                    if fields[-10].endswith('...'):
                        desc_end = -11 if len(fields) > 10 else -10
                    else:
                        desc_end = -10

                    desc = ' '.join(fields[:desc_end])

                    # Extract scientific name from first two words
                    desc_words = desc.split()
                    if len(desc_words) >= 2:
                        sci_name = desc_words[0] + ' ' + desc_words[1]
                    elif len(desc_words) == 1:
                        sci_name = desc_words[0]
                    else:
                        sci_name = "Unknown"

                    hit = {
                        'ReadName': current_read,
                        'Description': desc,
                        'ScientificName': sci_name,
                        'CommonName': common_name,
                        'Taxid': taxid,
                        'MaxScore': max_score,
                        'TotalScore': total_score,
                        'QueryCover': query_cover,
                        'Evalue': evalue,
                        'PercIdent': perc_ident,
                        'AccLen': acc_len,
                        'Accession': accession
                    }

                    hits_by_read[current_read].append(hit)
                    total_hits += 1

        # Check last query
        if current_read and current_read not in hits_by_read:
            hits_by_read[current_read] = []
            no_hit_queries += 1

    # Define header
    header = ['ReadName', 'Description', 'ScientificName', 'CommonName',
              'Taxid', 'MaxScore', 'TotalScore', 'QueryCover', 'Evalue',
              'PercIdent', 'AccLen', 'Accession']

    # Write all hits file
    with open(all_hits_file, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for read_name in sorted(hits_by_read.keys()):
            if not hits_by_read[read_name]:
                # No hits for this read - write NA entry
                na_hit = create_na_hit(read_name)
                f.write('\t'.join([na_hit[col] for col in header]) + '\n')
            else:
                for hit in hits_by_read[read_name]:
                    f.write('\t'.join([hit[col] for col in header]) + '\n')

    # Sort and write top 5 and top 1 for EACH READ
    with open(top5_file, 'w') as f5, open(top1_file, 'w') as f1:
        f5.write('\t'.join(header) + '\n')
        f1.write('\t'.join(header) + '\n')

        for read_name in sorted(hits_by_read.keys()):
            if not hits_by_read[read_name]:
                # No hits for this read - write NA entry to both files
                na_hit = create_na_hit(read_name)
                f5.write('\t'.join([na_hit[col] for col in header]) + '\n')
                f1.write('\t'.join([na_hit[col] for col in header]) + '\n')
            else:
                # Sort hits by score (highest first)
                sorted_hits = sorted(hits_by_read[read_name],
                                   key=score_hit,
                                   reverse=True)

                # Write top 5 for this read
                for hit in sorted_hits[:5]:
                    f5.write('\t'.join([hit[col] for col in header]) + '\n')

                # Write top 1 for this read
                if sorted_hits:
                    hit = sorted_hits[0]
                    f1.write('\t'.join([hit[col] for col in header]) + '\n')

    return all_hits_file, top5_file, top1_file, total_queries, total_hits, no_hit_queries

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python3 parse_blast_correct.py <input_file> <output_basename>")
        print("\nExample: python3 parse_blast_correct.py DEN0001.txt DEN0001")
        print("\nGenerates:")
        print("  - <basename>_all.tsv   (all hits for all reads)")
        print("  - <basename>_top5.tsv  (top 5 hits per read)")
        print("  - <basename>_top1.tsv  (best hit per read)")
        print("\nNote: Reads with no significant similarity get NA entries")
        sys.exit(1)

    input_file = sys.argv[1]
    output_base = sys.argv[2]

    all_file, top5_file, top1_file, num_reads, num_hits, no_hits = parse_blast_file(input_file, output_base)

    print(f"\n✓ Successfully parsed: {input_file}")
    print(f"\nSample statistics:")
    print(f"  - Total reads (queries): {num_reads}")
    print(f"  - Reads with hits: {num_reads - no_hits}")
    print(f"  - Reads with NO hits: {no_hits}")
    print(f"  - Total hits parsed: {num_hits}")
    if num_reads - no_hits > 0:
        print(f"  - Average hits per read (with hits): {num_hits/(num_reads - no_hits):.1f}")

    print(f"\nGenerated files:")
    print(f"  1. {all_file}")
    print(f"  2. {top5_file}")
    print(f"  3. {top1_file}")

    # Count rows in output files
    with open(all_file, 'r') as f:
        all_count = sum(1 for line in f) - 1
    with open(top5_file, 'r') as f:
        top5_count = sum(1 for line in f) - 1
    with open(top1_file, 'r') as f:
        top1_count = sum(1 for line in f) - 1

    print(f"\nOutput file rows (excluding header):")
    print(f"  - All hits: {all_count} rows")
    print(f"  - Top 5: {top5_count} rows")
    print(f"  - Top 1: {top1_count} rows (1 per read, {no_hits} with NA)")
