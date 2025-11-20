import os
import gzip
import collections
import csv
from datetime import datetime

# How to use it?
# Create the environment
# conda create -n fastq_tools python=3.12 biopython
# Active the environment
# conda activate fastq_tools 
# Change the input and export file.name below as you wish
# Then run this script in the directory where you also put you files/input directory
# python script_check_duplicated_indentifiers_from_concats.py
# Deactive the environment
# conda deactivate

# --- Configuration ---
# IMPORTANT: Update these paths to match your folder structure
INPUT_DIR = "fastq_compressed_files/"
REPORT_FILENAME = "fastq_duplication_report.tsv"

def analyze_fastq_duplicates(input_path):
    """
    Reads a single FASTQ file (gzipped or not), counts the total reads, 
    and identifies the number of unique and duplicate read headers.
    
    Args:
        input_path (str): The full path to the FASTQ file.
        
    Returns:
        tuple: (total_reads, unique_headers, duplicate_count) or (0, 0, 0) on error.
    """
    
    # Select the appropriate file opener function (handles .gz or plain files)
    is_gzipped = input_path.endswith('.gz')
    opener = gzip.open if is_gzipped else open
    
    read_id_counts = collections.defaultdict(int)
    total_reads = 0
    
    try:
        # Open file in text mode ('rt')
        with opener(input_path, 'rt') as infile: 
            
            while True:
                # 1: Header line (starts with '@')
                header_line = infile.readline()
                if not header_line:
                    break # End of file

                # Read and discard the next three lines (sequence, separator, quality scores)
                for _ in range(3):
                    if not infile.readline():
                        break 
                
                if not header_line.startswith('@'):
                    # Skip non-standard header lines if the file is corrupted/truncated
                    continue 

                total_reads += 1
                
                # --- Read ID Extraction Logic ---
                
                # Isolate the base read ID (everything after '@' up to the first space/tab/newline)
                end_of_id = header_line.find(' ')
                if end_of_id == -1:
                    end_of_id = len(header_line.strip())
                
                # Extract the base ID without '@'
                base_id = header_line[1:end_of_id]

                # Track the count for this base ID
                read_id_counts[base_id] += 1
        
    except Exception as e:
        print(f"ERROR: Failed to process {os.path.basename(input_path)}. Reason: {e}. Skipping file.")
        return 0, 0, 0
    
    # Calculate summary statistics
    unique_headers = len(read_id_counts)
    # The count of duplicated headers is the number of reads minus the number of unique IDs
    duplicate_count = total_reads - unique_headers 
    
    return total_reads, unique_headers, duplicate_count


def generate_duplication_report():
    """
    Main function to orchestrate the analysis and generate the report.
    """
    
    if not os.path.exists(INPUT_DIR):
        print(f"Error: Input directory '{INPUT_DIR}' not found. Please create it and place your FASTQ files inside.")
        return

    report_data = []
    
    # Iterate over all files in the input directory
    print(f"Starting analysis of files in: {INPUT_DIR}")
    
    for filename in os.listdir(INPUT_DIR):
        # Process only FASTQ files (gzipped or not)
        if filename.endswith(".fastq") or filename.endswith(".fastq.gz"):
            input_path = os.path.join(INPUT_DIR, filename)
            
            print(f"Analyzing {filename}...")
            total, unique, duplicate = analyze_fastq_duplicates(input_path)
            
            report_data.append({
                "Sample Name": filename,
                "Total Reads": total,
                "Unique Read Headers": unique,
                "Duplicate Read Headers": duplicate,
            })
            
    # Write the results to a TSV file
    if report_data:
        try:
            with open(REPORT_FILENAME, 'w', newline='') as tsvfile:
                fieldnames = ["Sample Name", "Total Reads", "Unique Read Headers", "Duplicate Read Headers"]
                writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')
                
                writer.writeheader()
                writer.writerows(report_data)
                
            print("\n------------------------------------------------------")
            print(f"Successfully generated duplication report: {REPORT_FILENAME}")
            print("------------------------------------------------------")
            
        except Exception as e:
            print(f"ERROR: Could not write report file. Reason: {e}")
    else:
        print("\nNo FASTQ files found to analyze.")

if __name__ == "__main__":
    generate_duplication_report()

