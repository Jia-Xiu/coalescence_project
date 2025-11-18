#!/bin/bash

# Set paths
input_csv="16S_mixing_exp._sample_list_30-04-2025.csv"
cleaned_dir="results_barbell_sg/results_with_adapters/cleaned_with_adapters_concats"
renamed_dir="results_barbell_sg/results_with_adapters/renamed_with_adapters_concats"
summary_dir="results_barbell_sg/results_with_adapters/summary_with_adapters_concats"
seqkit="/home/groups/VEO/tools/seqkit/v2.8.2/seqkit"

# Create output folders if they don’t exist
mkdir -p "$renamed_dir"
mkdir -p "$summary_dir"

# STEP 1: Rename samples

for i in `cat 16S_mixing_exp._sample_list_30-04-2025.csv`;do a=`echo $i | cut -d ',' -f 1`;b=`echo $i | cut -d ',' -f 2,3 | tr -d '\r' | sed "s/,/_/g"`; cp "${cleaned_dir}/${b}.fastq" "${renamed_dir}/${a}.fastq";echo $a;done

echo "✅ Renaming complete."

# STEP 2: Summarize read counts and length distributions
echo -e "\n📊 Summarizing reads and read lengths..."
summary_file="${summary_dir}/read_number_summary.tsv"
echo -e "Sample\tTotal_Reads" > "$summary_file"

for fastq in "$renamed_dir"/*.fastq; do
  sample=$(basename "$fastq" .fastq)
  total_reads=$($seqkit stats "$fastq" | awk 'NR==2 {print $4}')
  echo -e "${sample}\t${total_reads}" >> "$summary_file"

  # Generate read length distribution
  $seqkit fx2tab -l -n -i "$fastq" | awk '{print $2}' | sort | uniq -c | sort -nr > "${summary_dir}/length_distribution_${sample}.txt"
done

echo "✅ Summary written to $summary_file"



