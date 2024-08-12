#!/bin/bash
#SBATCH --job-name qiime2_sliva_classifier_training
#SBATCH --partition=gpu
#SBATCH --output tmp/qiime2_sliva_classifier_training.%j.out
#SBATCH --error  tmp/qiime2_sliva_classifier_training.%j.err
#SBATCH --mem=96G
#SBATCH --cpus-per-task=12
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=xiu.jia@uni-jena.de



#-- Command section ------------------------

source /vast/groups/VEO/tools/miniconda3_2024/etc/profile.d/conda.sh && conda activate qiime2-amplicon-2024.5

mkdir -p silva_db
cd silva_db

# follow this tutorial https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494


# https://www.arb-silva.de/projects/ssu-ref-nr/
# downloading the relevant taxonomy and sequence files from the SILVA.
# Import these files into QIIME 2.
qiime rescript get-silva-data \
    --p-version '138.1' \
    --p-target 'SSURef_NR99' \
    --o-silva-sequences silva-138.1-ssu-nr99-rna-seqs.qza \
    --o-silva-taxonomy silva-138.1-ssu-nr99-tax.qza

# Import these files into QIIME 2.
qiime rescript reverse-transcribe \
    --i-rna-sequences silva-138.1-ssu-nr99-rna-seqs.qza \
    --o-dna-sequences silva-138.1-ssu-nr99-seqs.qza

# “Culling” low-quality sequences with cull-seqs
# Here we’ll remove sequences that contain 5 or more ambiguous bases (IUPAC compliant ambiguity bases) and any homopolymers that are 8 or more bases in length. These are the default parameters. 
qiime rescript cull-seqs \
    --i-sequences silva-138.1-ssu-nr99-seqs.qza \
    --o-clean-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza

# Filtering sequences by length and taxonomy
# remove rRNA gene sequences that do not meet the following criteria: Archaea (16S) >= 900 bp, Bacteria (16S) >= 1200 bp, and any Eukaryota (18S) >= 1400 bp.
qiime rescript filter-seqs-length-by-taxon \
    --i-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza \
    --i-taxonomy silva-138.1-ssu-nr99-tax.qza \
    --p-labels Archaea Bacteria Eukaryota \
    --p-min-lens 900 1200 1400 \
    --o-filtered-seqs silva-138.1-ssu-nr99-seqs-filt.qza \
    --o-discarded-seqs silva-138.1-ssu-nr99-seqs-discard.qza

# Dereplication of sequences and taxonomy
qiime rescript dereplicate \
    --i-sequences silva-138.1-ssu-nr99-seqs-filt.qza  \
    --i-taxa silva-138.1-ssu-nr99-tax.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
    --o-dereplicated-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza


# Make amplicon-region specific classifier
qiime feature-classifier extract-reads \
    --i-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza \
    --p-f-primer AGRGTTYGATYMTGGCTCAG \
    --p-r-primer CGGYTACCTTGTTACGACTT \
    --p-n-jobs 2 \
    --p-read-orientation 'forward' \
    --o-reads silva-138.1-ssu-nr99-seqs-27f-1492r.qza

# Dereplication of sequences and taxonomy
qiime rescript dereplicate \
    --i-sequences silva-138.1-ssu-nr99-seqs-27f-1492r.qza \
    --i-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza \
    --p-mode 'uniq' \
    --o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-27f-1492r-uniq.qza \
    --o-dereplicated-taxa  silva-138.1-ssu-nr99-tax-27f-1492r-derep-uniq.qza

# Make amplicon-region specific classifier
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads silva-138.1-ssu-nr99-seqs-27f-1492r-uniq.qza \
    --i-reference-taxonomy silva-138.1-ssu-nr99-tax-27f-1492r-derep-uniq.qza \
    --o-classifier silva-138.1-ssu-nr99-27f-1492r-classifier.qza


conda deactivate


