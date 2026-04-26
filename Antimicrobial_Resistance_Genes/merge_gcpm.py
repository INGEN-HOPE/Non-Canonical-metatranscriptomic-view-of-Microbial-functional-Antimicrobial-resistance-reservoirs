import glob
import os

def merge_files():
    # Get all gcpm files, excluding summary statistics
    files = glob.glob("*_gcpm.csv")
    files = [f for f in files if "summary_statistics" not in f]
    
    # Create output file and write header from first file
    with open(files[0], 'r') as first_file:
        header = first_file.readline().strip()
    
    with open("merged_gcpm.csv", 'w') as outfile:
        # Write header
        outfile.write(f"Sample,{header}\n")
        
        # Process each file
        for file in files:
            sample_name = file.replace("_gcpm.csv", "")
            
            # Skip header and add sample name to each line
            with open(file, 'r') as infile:
                next(infile)  # skip header
                for line in infile:
                    outfile.write(f"{sample_name},{line}")

if __name__ == "__main__":
    merge_files()
    print("Files merged successfully into merged_gcpm.csv")
