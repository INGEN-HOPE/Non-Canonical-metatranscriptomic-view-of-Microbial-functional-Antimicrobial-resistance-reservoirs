import pandas as pd
import glob
import os

def combine_rgi_outputs(input_dir, output_prefix):
    """
    Combine RGI CARD output files from multiple samples into consolidated reports.
    
    Args:
        input_dir (str): Directory containing RGI output files
        output_prefix (str): Prefix for output files
    """
    # Initialize empty lists to store dataframes
    gene_mapping_dfs = []
    allele_mapping_dfs = []
    
    # Get all gene and allele mapping files
    gene_files = glob.glob(os.path.join(input_dir, "*.gene_mapping_data_filtered.txt"))
    allele_files = glob.glob(os.path.join(input_dir, "*.allele_mapping_data_filtered.txt"))
    
    # Process gene mapping files
    for file in gene_files:
        sample_name = os.path.basename(file).split("_rgi_output")[0]
        try:
            df = pd.read_csv(file, sep='\t')
            df['Sample_ID'] = sample_name
            gene_mapping_dfs.append(df)
        except Exception as e:
            print(f"Error processing gene mapping file {file}: {str(e)}")
    
    # Process allele mapping files
    for file in allele_files:
        sample_name = os.path.basename(file).split("_rgi_output")[0]
        try:
            df = pd.read_csv(file, sep='\t')
            df['Sample_ID'] = sample_name
            allele_mapping_dfs.append(df)
        except Exception as e:
            print(f"Error processing allele mapping file {file}: {str(e)}")
    
    # Combine all dataframes
    if gene_mapping_dfs:
        combined_gene_df = pd.concat(gene_mapping_dfs, ignore_index=True)
        # Move Sample_ID to first column
        cols = ['Sample_ID'] + [col for col in combined_gene_df.columns if col != 'Sample_ID']
        combined_gene_df = combined_gene_df[cols]
        # Save combined gene mapping data
        output_gene_file = f"{output_prefix}_combined_gene_mapping.txt"
        combined_gene_df.to_csv(output_gene_file, sep='\t', index=False)
        print(f"Combined gene mapping data saved to: {output_gene_file}")
    
    if allele_mapping_dfs:
        combined_allele_df = pd.concat(allele_mapping_dfs, ignore_index=True)
        # Move Sample_ID to first column
        cols = ['Sample_ID'] + [col for col in combined_allele_df.columns if col != 'Sample_ID']
        combined_allele_df = combined_allele_df[cols]
        # Save combined allele mapping data
        output_allele_file = f"{output_prefix}_combined_allele_mapping.txt"
        combined_allele_df.to_csv(output_allele_file, sep='\t', index=False)
        print(f"Combined allele mapping data saved to: {output_allele_file}")

    return {
        'gene_mapping': f"{output_prefix}_combined_gene_mapping.txt" if gene_mapping_dfs else None,
        'allele_mapping': f"{output_prefix}_combined_allele_mapping.txt" if allele_mapping_dfs else None
    }

if __name__ == "__main__":
    # Example usage
    input_directory = "."  # Current directory
    output_prefix = "combined_rgi"
    output_files = combine_rgi_outputs(input_directory, output_prefix)
