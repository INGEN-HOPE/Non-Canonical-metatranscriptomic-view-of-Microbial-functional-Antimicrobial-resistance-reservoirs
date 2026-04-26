import pandas as pd
import numpy as np
import glob
import os

def calculate_gcpm(allele_mapping_file):
    """
    Calculate Gene Coverage Per Million (GCPM) values from RGI allele mapping output
    using the formula: (counts/gene length) × 10^6 / Σ(counts/gene length) for all n genes

    Parameters:
    allele_mapping_file: Path to the allele_mapping_data.txt from RGI output

    Returns:
    DataFrame with GCPM values for each ARG
    """
    # Read allele mapping data
    allele_data = pd.read_csv(allele_mapping_file, sep='\t')

    # Extract relevant columns
    gene_data = allele_data[[
        'Reference Sequence',
        'ARO Term',
        'Reference Length',
        'Completely Mapped Reads',
        'AMR Gene Family',
        'Drug Class',
        'Resistance Mechanism',
        'Resistomes & Variants: Observed Pathogen(s)',
        'Percent Coverage'
    ]].copy()

    # Calculate counts/gene length for each gene
    gene_data['counts_per_length'] = gene_data['Completely Mapped Reads'] / gene_data['Reference Length']

    # Calculate the sum of counts_per_length for normalization
    total_counts_per_length = gene_data['Completely Mapped Reads'].sum()

    # Calculate GCPM
    if total_counts_per_length > 0:  # Avoid division by zero
        gene_data['gcpm'] = (gene_data['counts_per_length'] * 1e6) / (total_counts_per_length / gene_data['Reference Length']) 
    else:
        gene_data['gcpm'] = 0

    # Clean up intermediate calculation column
    gene_data = gene_data.drop('counts_per_length', axis=1)

    # Set index to Reference Sequence
    gene_data = gene_data.set_index('Reference Sequence')

    return gene_data

def process_multiple_samples(directory=".", output_dir="gcpm"):
    """
    Process multiple samples and create individual GCPM files for each sample

    Parameters:
    directory: Directory containing RGI output files
    output_dir: Directory where GCPM files will be saved

    Returns:
    Dictionary with summary statistics for each sample
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Find all allele mapping files
    allele_files = glob.glob(os.path.join(directory, "*.allele_mapping_data_filtered.txt"))

    summary_stats = {}

    for allele_file in allele_files:
        # Extract sample name from file name
        sample_name = os.path.basename(allele_file).replace(".allele_mapping_data_filtered.txt", "")

        print(f"Processing sample: {sample_name}")

        try:
            # Calculate GCPM for this sample
            gcpm_values = calculate_gcpm(allele_file)

            # Save individual sample GCPM file
            output_file = os.path.join(output_dir, f"{sample_name}_gcpm.csv")
            gcpm_values.to_csv(output_file)

            # Calculate summary statistics
            gcpm_series = gcpm_values['gcpm']
            summary_stats[sample_name] = {
                'total_args': (gcpm_series > 0).sum(),
                'max_gcpm': gcpm_series.max(),
                'mean_gcpm': gcpm_series.mean()
            }

            print(f"GCPM values saved to {output_file}")

        except Exception as e:
            print(f"Error processing {sample_name}: {str(e)}")
            continue

    # Create summary file
    summary_df = pd.DataFrame.from_dict(summary_stats, orient='index')
    summary_df.to_csv(os.path.join(output_dir, "summary_statistics.csv"))

    return summary_stats

if __name__ == "__main__":
    try:
        # Process all samples and save individual GCPM files
        print("Starting GCPM calculation for all samples...")
        summary_stats = process_multiple_samples()

        # Print summary statistics
        print("\nSummary of GCPM values:")
        for sample, stats in summary_stats.items():
            print(f"\nSample: {sample}")
            print(f"Total ARGs detected: {stats['total_args']}")
            print(f"Max GCPM value: {stats['max_gcpm']:.2f}")
            print(f"Mean GCPM value: {stats['mean_gcpm']:.2f}")

        print(f"\nAll GCPM files have been saved in the 'gcpm' directory")
        print("Summary statistics have been saved to 'gcpm/summary_statistics.csv'")

    except Exception as e:
        print(f"An error occurred: {str(e)}")
