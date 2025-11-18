# Pipeline to analyze raw reads of full length 16S RNA gene
We got **14.29 M raw reads** of full length 16S RNA gene in pod5 files, which were generated from the FLO-MIN114 flowcell :dna: using the SQK-LSK114 kit.

## 1. Basecalling
I used [**Dorado (v0.9.1)**](https://github.com/nanoporetech/dorado/tree/release-v0.9) to do basecalling.
We got 17M reads after basecalling by Dorado, with 4M “redundant” reads, which are so called "simplex have duplex offsprings". I used non_redundant.fastq.gz (12.9M reads) file for further demultiplexing test. \
Xiu: Maybe I should do quality and sequence length control before demultiplexing.

```
/.../dorado/v0.9.1/bin/dorado  basecaller --emit-fastq \
        /.../dorado/v0.9.1/models/dna_r10.4.1_e8.2_400bps_sup@v4.3.0 \
        --no-trim \
        ../raw_data_jena/240202_16S_amplicons_Xiu/20240202_1646_MN41792_FAW73518_3914cd2c/pod5/ \
        -o results_dorado_0.9_with_adpters

## nanoplot
source /vast/groups/VEO/tools/anaconda3/etc/profile.d/conda.sh && conda activate nanoplot_v1.41.3

NanoPlot -t 2 --fastq results_dorado_0.9_with_adpters/calls_2025-03-28_T19-44-13.fastq -o results_nanoplot/Dorado_0.9_with_adapter/all

/.../dorado/v0.9.1/bin/dorado  basecaller --emit-fastq \
        /.../dorado/v0.9.1/models/dna_r10.4.1_e8.2_400bps_sup@v4.3.0 \
        --trim 'adapters' \
        ../raw_data_jena/240202_16S_amplicons_***/20240202_1646_MN41792_FAW73518_3914cd2c/pod5/ \
        -o results_dorado_0.9

gzip -c results_dorado_0.9/calls_2025-02-25_T15-31-26.fastq.gz > results_dorado_0.9/calls_2025-02-25_T15-31-26.fastq


# make nanoplot for all the samples
NanoPlot -t 2 --fastq results_dorado_0.9/calls_2025-02-25_T15-31-26.fastq -o results_nanoplot/Dorado_0.9/all

```
NanoStats
```
General summary:         
Mean read length:                  1,535.2
Mean read quality:                     9.6
Median read length:                1,549.0
Median read quality:                  14.3
Number of reads:              14,565,365.0
Read length N50:                   1,555.0
STDEV read length:                 1,373.3
Total bases:              22,361,283,770.0
```

## 2. Demultiplexing
### 2.1. Barbell
We used [**barbell-sg-v0.1.5**](https://github.com/rickbeeloo/barbell) which is developed by **Rick Beloo** to demultiplex raw 16S reads. \
From 14.56 M raw reads, we got **6,267,380 reads** assigned family taxonomic level.
```
#!/bin/bash
#SBATCH --job-name demultiplex_barbell-sg-v0.1.5_with_adapters
#SBATCH --partition=short
#SBATCH --mem=4G
#SBATCH --cpus-per-task=8
#SBATCH --output tmp/demultiplex_barbell-sg-v0.1.5_with_adapters.%j.out
#SBATCH --error  tmp/demultiplex_barbell-sg-v0.1.5_with_adapters.%j.err
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=***


# Run barbell-sg using the following path

barbell="/.../my_tools/barbell-sg-v0.1.5/target/release/barbell"


# Define paths of main output folder and input files
out_dir="results_barbell_sg/results_with_adapters"
raw_fastq="results_dorado_0.9_with_adpters/calls_2025-03-28_T19-44-13.fastq"
fwd_primers="results_barbell_sg/forward_primer_barcodes.fasta"
rev_primers="results_barbell_sg/reverse_primer_barcodes.fasta"

# Create output directory if it doesn't exist
mkdir -p "$out_dir"

# Annotate
$barbell annotate \
    -i "$raw_fastq" \
    -q "$fwd_primers","$rev_primers" \
    -o $out_dir/annotations_with_adapters.txt \
    -t 8 \
    --tune

# Inspect
$barbell inspect -i $out_dir/annotations_with_adapters.txt

# Option 1: consider reads with only 2 barcode annotations

# Filter
$barbell filter \
    -i $out_dir/annotations_with_adapters.txt \
    -o $out_dir/filtered_with_adapters_2barcodes.txt \
    -f results_barbell_sg/rapid_filters_250_1300_1500.txt

# Trim
$barbell trimm \
    -i $out_dir/filtered_with_adapters_2barcodes.txt \
    -r "$raw_fastq" \
    -o $out_dir/trimmed_with_adapters_2barcodes

# Option 2: consider concatenated reads

# Filter
$barbell filter \
    -i $out_dir/annotations_with_adapters.txt \
    -o $out_dir/filtered_with_adapters_conctas.txt \
    -f results_barbell_sg/rapid_filters_250_1300_1500_2_3_concats.txt

# Trim
$barbell trimm \
    -i $out_dir/filtered_with_adapters_conctas.txt \
    -r "$raw_fastq" \
    -o $out_dir/trimmed_with_adapters_conctas
```

### 2.2. Rename demultiplexed fastq files


## 3. Taxonomic assignment

### 3.1 EMU
We assigned taxonomy to the raw reads by using [**Emu**: species-level taxonomic abundance for full-length 16S reads](https://github.com/treangenlab/emu).

* Curry, K.D., Wang, Q., Nute, M.G. et al. Emu: species-level microbial community profiling of full-length 16S rRNA Oxford Nanopore sequencing data. Nat Methods 19, 845–853 (2022). https://doi.org/10.1038/s41592-022-01520-4

run this code: `bash script_06_taxonomy_annotation_emu.sh`:
```
#!/bin/bash
mkdir -p tmp/emu_taxonomy_annotation emu_out

for fastq in results_barbell_sg/results_with_adapters/renamed_with_adapters_concats/*.fastq; do
    sample=$(basename "$fastq" .fastq)
    sbatch --job-name="$sample" sbatch_06_taxonomy_annotation_emu.sh "$fastq"
done
```

sbatch_06_taxonomy_annotation_emu.sh:
```
#!/bin/bash
#SBATCH --job-name=emu_anno_%x
#SBATCH --partition=standard
#SBATCH --output=tmp/emu_taxonomy_annotation/%x.%j.out
#SBATCH --error=tmp/emu_taxonomy_annotation/%x.%j.err
#SBATCH --mem=10G
#SBATCH --cpus-per-task=10

# Load emu in conda environment
source /vast/groups/VEO/tools/miniconda3_2024/etc/profile.d/conda.sh 
conda activate emu_v3.4.5 
emu -h #database /veodata/03/databases/emu/v240807

# Input file is passed as argument
fastq_file_name=$1

# Create output directory (unique per file)
sample_name=$(basename "$fastq_file_name" .fastq)
output_dir="emu_out"
mkdir -p "$output_dir"

emu abundance "$fastq_file_name" \
        --db /veodata/03/databases/emu/v240807 \
        --keep-counts \
        --keep-read-assignments \
        --threads 10 \
        --output-dir  "$output_dir"

emu combine-outputs emu_out/ species       

exit
```

### 3.2 Kraken2
We assigned taxonomy to the raw reads by using [**Kraken2**](https://github.com/DerrickWood/kraken2/wiki/Manual). To report the output, I used Braken, a [customized python script](https://combine_kreports.py)  by jennifer.lu717@gmail.com. See here: https://github.com/jenniferlu717/Bracken?tab=readme-ov-file and https://ccb.jhu.edu/software/bracken/
```
#!/bin/bash
#SBATCH --job-name kraken2_barbell_raw_reads
#SBATCH --partition=standard
#SBATCH --output tmp/kraken2_barbell_raw_reads.%j.out
#SBATCH --error  tmp/kraken2_barbell_raw_reads.%j.err
#SBATCH --mem=20G
#SBATCH --cpus-per-task=20
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=...

# https://github.com/DerrickWood/kraken2/wiki/Manual
kraken2_db="/.../database/16S_SILVA138_k2db/"

# Input directory containing the FASTQ files
input_dir="results_barbell_sg/results_with_adapters/renamed_with_adapters_concats/"

# Output directory where filtered FASTQ files will be saved
output_dir="results_kraken2_barbell_adapter_concats_silva_182"


# Assign taxonomy

mkdir -p "$output_dir"

for file in $input_dir/*.fastq; do
  # Extract the base name of the file (without path and extension)

  filename=$(basename "$file" .fastq)

  # Run Kraken2 on the current FASTQ file
  /.../tools/kraken2/v2.1.5/kraken2 \
          --db $kraken2_db \
          --threads 10 \
          --use-mpa-style \
          --report $output_dir/${filename}_report.txt \
          --output $output_dir/${filename}_output.txt $file

  echo "Processed $file with Kraken2 and saved results to $output_dir"

done


echo "\nstart Kraken2 result summary by python code"

# activate the python env

source /.../my_tools/mypyenv/bin/activate

python3.9 /.../my_tools/KrakenTools/combine_kreports.py -r $output_dir/*_report.txt -o kraken2_combined_abundance_report_adapter_concats_silva_182.txt

deactivate
```

Now enjoy the downstream analysis :sparkles:

### 3.3 Other analysis attempt
I also tried **QIIME2** platform for sequence analysis.\
By using VSEARCH from QIIME2, I found majority reads are unique because the higher error rates of the reads. From 3.6M reads, only 4263 reads can be rereplicated (sea the table from the table.qzv file).\
Current available denoising approaches, such as DADA2, are inappropriate for Nanopore long reads. Because current ONT data has too high an error rate for the DADA2 approach to be valid. DADA2 currently supports PacBio circular consensus sequencing but not nanopore reads (Callahan et al. 2019). Ref: https://github.com/benjjneb/dada2/issues/759 and https://github.com/benjjneb/dada2/issues/1364 

## 4. Prepare reads for submission
As we also cut and kept concatenated reads from Barbell for further analysis, some reads in FASTQs have duplicated indentifier. These reads cannot be accepted by NCBI or ENA. So I added extension to these duplicated reads.

```
import os
import gzip
import sys

# --- Configuration ---
# IMPORTANT: Update these paths to match your folder structure on your machine
INPUT_DIR = "renamed_with_adapters_concats/"
OUTPUT_DIR = "renamed_with_adapters_concats_header_updated/"


def rename_duplicate_fastq_headers(input_path, output_path):
    """
    Reads a single FASTQ (or gzipped FASTQ) file, renames only the duplicate 
    read headers using the _DUP_N convention, and writes the output.
    """
    
    # Select the appropriate file opener function (handles .gz or plain files)
    is_gzipped = input_path.endswith('.gz')
    opener = gzip.open if is_gzipped else open
    
    read_id_counts = {}
    reads_renamed = 0
    total_reads = 0
    
    try:
        # Open files using the determined opener function
        with opener(input_path, 'rt') as infile, \
             opener(output_path, 'wt') as outfile: 
            
            while True:
                # 1: Header line (starts with '@')
                header_line = infile.readline()
                if not header_line:
                    break # End of file

                # 2-4: Read the rest of the FASTQ block
                sequence_line = infile.readline()
                quality_sep_line = infile.readline()
                quality_line = infile.readline()
                
                if not (sequence_line and quality_sep_line and quality_line):
                    break 

                total_reads += 1
                
                # --- Read ID Extraction and Renaming Logic ---
                
                # Isolate the base read ID (everything after '@' up to the first space/tab/newline)
                end_of_id = header_line.find(' ')
                if end_of_id == -1:
                    end_of_id = len(header_line.strip())
                
                original_id_full = header_line[:end_of_id] # e.g., @read_id
                base_id = original_id_full[1:]              # e.g., read_id

                # Track the count for this base ID
                read_id_counts[base_id] = read_id_counts.get(base_id, 0) + 1
                current_count = read_id_counts[base_id]
                
                new_header_line = header_line
                
                if current_count > 1:
                    # This is a duplicate. Rename it.
                    reads_renamed += 1
                    
                    # Calculate the duplicate index (1 for the first duplicate, 2 for the second, etc.)
                    duplicate_index = current_count - 1 
                    
                    # Construct the new unique ID (e.g., ID_DUP_1)
                    suffix = f"_DUP_{duplicate_index}"
                    new_base_id_full = original_id_full + suffix
                    
                    # Replace the old base ID with the new one in the header line
                    new_header_line = header_line.replace(original_id_full, new_base_id_full, 1)

                # Write the modified block
                outfile.write(new_header_line)
                outfile.write(sequence_line)
                outfile.write(quality_sep_line)
                outfile.write(quality_line)
        
    except Exception as e:
        print(f"ERROR: Failed to process {input_path}. Reason: {e}")
        if os.path.exists(output_path):
            os.remove(output_path)
        return

    print(f"-> FINISHED: {os.path.basename(input_path)} | Total reads: {total_reads} | Renamed: {reads_renamed}")


# --- Main Execution Loop ---

# 1. Ensure the output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)
print(f"Created/verified output directory: {OUTPUT_DIR}\n")

# 2. Iterate over all files in the input directory
print(f"Starting batch processing of files in: {INPUT_DIR}")
for filename in os.listdir(INPUT_DIR):
    # Process only FASTQ files (gzipped or not)
    if filename.endswith(".fastq") or filename.endswith(".fastq.gz"):
        input_path = os.path.join(INPUT_DIR, filename)
        output_path = os.path.join(OUTPUT_DIR, filename)
        
        print(f"Processing {filename}...")
        rename_duplicate_fastq_headers(input_path, output_path)

print("\nBatch processing complete. The corrected files are in the 'files_header_updated/' folder.")
```
Check one file that have been re-named successfully
`cat TEcontrol_r2.fastq | awk '/^@.*_DUP_/{for(i=0;i<4;i++) {print; getline}}`

