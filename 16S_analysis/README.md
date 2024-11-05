# Pipeline to analyze raw reads of full length 16S RNA gene
We got 14.29 M reads of full length 16S RNA gene in pod5 files, which were generated from the FLO-MIN114 flowcell :dna: using the SQK-LSK114 kit.

## Basecalling
I used [Dorado (v0.6.0)](https://github.com/nanoporetech/dorado/tree/release-v0.6.0) to do basecalling.
We got 17M reads after basecalling by Dorado, with 4M “redundant” reads, which are so called "simplex have duplex offsprings". I used non_redundant.fastq.gz (12.9M reads) file for further demultiplexing test. \
Xiu: Maybe I should do quality and sequence length control before demultiplexing.
```
/.../dorado/v0.6.0/bin/dorado duplex \
/.../dorado/v0.6.0/models/dna_r10.4.1_e8.2_400bps_sup@v4.3.0 \
data/pod5/ > duplex.bam

/.../tools/samtools/v1.17/bin/samtools view --tag dx:1 --tag dx:0 \
results_dorado_no_demultiplexing/duplex.bam | gzip -9 > results_dorado_no_demultiplexing/non_redundant.fastq.gz
```
## Demultiplexing
### Barbell
We used [Barbell](https://github.com/rickbeeloo/barbell) which is developed by Rick Beloo to demultiplex raw 16S reads. \
From 12.9M raw reads, we got 3,686,727 reads assigned family taxonomic level.
```
# active the env
source /.../anaconda3/etc/profile.d/conda.sh
conda activate Barbell_new

# run barbell
/.../my_tools/barbell-main_2/target/release/barbell \
        -c configs_barbell/config.toml \
        -s configs_barbell/samples.txt \
        -r results_dorado/non_redundant.fastq \
        -t 16 \
        -o results_barbell
```

### prob-edit-rs
I will also use [prob-edit-rs](https://github.com/rickbeeloo/prob-edit-rs) for demultiplexing as it might be more precise than Barbell.

## Taxonomic assignment
### Kraken2
We assigned taxonomy to the raw reads by using [Kraken2](https://github.com/DerrickWood/kraken2/wiki/Manual). To report the output, I used Braken, a [customized python script](https://combine_kreports.py)  by jennifer.lu717@gmail.com. See here: https://github.com/jenniferlu717/Bracken?tab=readme-ov-file and https://ccb.jhu.edu/software/bracken/
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
kraken2_db="/work/groups/VEO/databases/kraken2/v20180901"

# Input directory containing the FASTQ files
input_dir="results_barbell"

# Output directory where filtered FASTQ files will be saved
output_dir="results_kraken2_barbell_2"


# Assign taxonomy

mkdir -p "$output_dir"

for file in $input_dir/*.fastq; do
  # Extract the base name of the file (without path and extension)

  filename=$(basename "$file")

  # Run Kraken2 on the current FASTQ file
  /.../VEO/tools/kraken2/v2.1.2/kraken2 \
          --db $kraken2_db \
          --threads 20 \
          --report $output_dir/${filename}_report.txt \
          --output $output_dir/${filename}_output.txt $file

  echo "Processed $file with Kraken2 and saved results to $kraken_output_dir"

done


echo "\nstart Kraken2 result summary by python code"

# activate the python env

source /.../my_tools/mypyenv/bin/activate

python3.9 /.../my_tools/KrakenTools/combine_kreports.py -r $output_dir/*_report.txt -o kraken2_combined_abundance_report.txt

deactivate
```

Now enjoy the downstream analysis :sparkles:


