# Pipeline: nanopore sequencing raw reads of full length 16S RNA gene
## Basecalling
I used Dorado (v0.6.0) to do base calling.
We got 17M reads after basecalling by Dorado, with 4M “redundant” simplex have duplex offsprings. I used non_redundant.fastq.gz (12M reads) file for further demultiplexing test. 
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
We used [Barbell](https://github.com/rickbeeloo/barbell) which is developed by Rick Beloo to demultiplex raw 16S reads
From 12,9M raw reads, we got 3,686,727 reads assigned family taxonomic level.
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
