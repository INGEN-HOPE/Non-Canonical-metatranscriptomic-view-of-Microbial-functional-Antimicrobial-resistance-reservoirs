import pandas as pd
import glob
import os

def filter_antiseptic_files(directory_pattern="*.allele_mapping_data.txt"):
    # Get all files matching the pattern
    files = glob.glob(directory_pattern)
    
    for file in files:
        try:
            # Read the data
            df = pd.read_csv(file, sep='\t')
            
            # Filter out rows where Drug Class is ONLY "disinfecting agents and antiseptics"
            filtered_df = df[df['Drug Class'] != 'disinfecting agents and antiseptics']
            
            # Create output filename
            output_file = file.replace('.txt', '_filtered.txt')
            
            # Write filtered data
            filtered_df.to_csv(output_file, sep='\t', index=False)
            
            print(f"Processed {file} -> {output_file}")
            print(f"Removed {len(df) - len(filtered_df)} rows")
            
        except Exception as e:
            print(f"Error processing {file}: {str(e)}")

# Run the function
filter_antiseptic_files()
