# ITS1 fungal metabarcoding pipeline (QIIME 2)
# -------------------------------------------------------------------
# Project: ITS1 – 5 island comparison
#
# QIIME 2: amplicon-2024.5 
# Reference database: UNITE v 9.0 
#
# ITS1-specific recommendations:
#   - Pauvert et al. 2019 "Bioinformatics matters" (ITS1, single-end R1)
#   - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8765055/
#
# Strategy:
#   - Only R1 reads (single-end; ITS1 highly variable in length)
#   - Primer trimming with cutadapt
#   - Denoising with DADA2 denoise-single
#
# Primers (EMP ITS1 protocol):
#   - ITS1F:  CTTGGTCATTTAGAGGAAGTAA
#   - ITS2:   GCTGCGTTCTTCATCGATGC (reverse, reverse-complement used for adapter)
#
# Project structure:
#   - Two Illumina runs: run1/ and run2/
#
# NOTE: This script assumes:
#   - raw_data/ directory exists inside each runN/ directory
#   - Conda environments:
#       qiime2-amplicon-2024.5
#       base (with seqkit installed)
# -------------------------------------------------------------------

############# ACTIVATE QIIME2 ENVIRONMENT ###########################

cd ~/qiime2/analysis/ITS_SRE/analysis_ITS

conda activate qiime2-amplicon-2024.5

mkdir -p run1
mkdir -p run2

############# DENOISING – RUN 1 #####################################

cd run1

# 1. Import single-end R1 reads
# -----------------------------

# Create folder for R1 reads only
mkdir -p raw_data_R1

# Copy only R1 FASTQ files from raw_data
cp ./raw_data/*_R1_* ./raw_data_R1/

# Import into QIIME 2 (Casava 1.8 single-end format)
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path raw_data_R1 \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-single-end.qza

# Demultiplexing summary
qiime demux summarize \
  --i-data demux-single-end.qza \
  --o-visualization demux-single-end.qzv

# qiime tools view demux-single-end.qzv

# 2. Cutadapt primer trimming
# ---------------------------
# We trim both the forward primer (ITS1F) and the reverse-complement of ITS2.

qiime cutadapt trim-single \
  --i-demultiplexed-sequences demux-single-end.qza \
  --p-cores 25 \
  --p-front CTTGGTCATTTAGAGGAAGTAA \
  --p-adapter GCATCGATGAAGAACGCAGC \
  --p-error-rate 0.1 \
  --p-overlap 3 \
  --p-match-adapter-wildcards True \
  --o-trimmed-sequences trimmed.qza

# Summarize trimmed reads
qiime demux summarize \
  --i-data trimmed.qza \
  --o-visualization trimmed.qzv

# qiime tools view trimmed.qzv

#####################################################
### Optional: check primer trimming with seqkit   ###
#####################################################

# Export sequences before and after cutadapt
qiime tools export --input-path demux-single-end.qza --output-path before_cutadapt
qiime tools export --input-path trimmed.qza           --output-path after_cutadapt

# Switch to base environment with seqkit
conda activate base

# Stats and primer counts with seqkit
# See: https://bioinf.shenwei.me/seqkit/usage/#grep

## Stats – R1 before cutadapt

# Combine all FASTQ files into a temporary file
cat before_cutadapt/*.fastq.gz > before_cutadapt/combined.fastq.gz

# Get summary statistics
seqkit stats before_cutadapt/combined.fastq.gz >> stats_primers.txt

# Remove temporary combined file
rm before_cutadapt/combined.fastq.gz

## Stats – R1 after cutadapt

cat after_cutadapt/*.fastq.gz > after_cutadapt/combined.fastq.gz
seqkit stats after_cutadapt/combined.fastq.gz >> stats_primers.txt
rm after_cutadapt/combined.fastq.gz

## Count primer occurrences (before and after trimming)

# Forward primer (ITS1F)
cat before_cutadapt/*.fastq.gz | \
  seqkit grep -s -i -d -p CTTGGTCATTTAGAGGAAGTAA -C >> count_primers.txt

cat after_cutadapt/*.fastq.gz | \
  seqkit grep -s -i -d -p CTTGGTCATTTAGAGGAAGTAA -C >> count_primers.txt

# Reverse primer (ITS2, reverse-complement)
cat before_cutadapt/*.fastq.gz | \
  seqkit grep -s -i -d -p GCATCGATGAAGAACGCAGC -C >> count_primers.txt

cat after_cutadapt/*.fastq.gz | \
  seqkit grep -s -i -d -p GCATCGATGAAGAACGCAGC -C >> count_primers.txt


# 3. Denoising / ASV inference (DADA2)
# ------------------------------------
# For ITS1 (highly variable length), we:
#   - run denoise-single
#   - set trunc-len to 0 (no truncation by position)
#   - use --p-max-ee 8.0 and --p-trunc-q 8 as recommended for ITS1

qiime dada2 denoise-single \
  --i-demultiplexed-seqs trimmed.qza \
  --p-trunc-len 0 \
  --p-trim-left 0 \
  --p-trunc-q 8 \
  --p-max-ee 8.0 \
  --p-n-threads 18 \
  --p-pooling-method pseudo \
  --output-dir dada2 \
  --verbose

#   --p-pooling-method pseudo increases sensitivity to shared ASVs

# Export DADA2 denoising statistics as TSV
qiime tools export \
  --input-path dada2/denoising_stats.qza \
  --output-path dada2

# Visualize DADA2 stats
qiime metadata tabulate \
  --m-input-file dada2/denoising_stats.qza \
  --o-visualization dada2/stats-dada2.qzv

# qiime tools view dada2/stats-dada2.qzv

# Summarize feature table
qiime feature-table summarize \
  --i-table dada2/table.qza \
  --o-visualization dada2/table.qzv

# qiime tools view dada2/table.qzv

cd ..

############# DENOISING – RUN 2 #####################################

cd run2

# 1. Import single-end R1 reads
# -----------------------------

mkdir -p raw_data_R1
cp ./raw_data/*_R1_* ./raw_data_R1/

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path raw_data_R1 \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-single-end.qza

qiime demux summarize \
  --i-data demux-single-end.qza \
  --o-visualization demux-single-end.qzv

# qiime tools view demux-single-end.qzv

# 2. Cutadapt primer trimming
# ---------------------------

qiime cutadapt trim-single \
  --i-demultiplexed-sequences demux-single-end.qza \
  --p-cores 25 \
  --p-front CTTGGTCATTTAGAGGAAGTAA \
  --p-adapter GCATCGATGAAGAACGCAGC \
  --p-error-rate 0.1 \
  --p-overlap 3 \
  --p-match-adapter-wildcards True \
  --o-trimmed-sequences trimmed.qza

qiime demux summarize \
  --i-data trimmed.qza \
  --o-visualization trimmed.qzv

# qiime tools view trimmed.qzv

#####################################################
### Optional: check primer trimming with seqkit   ###
#####################################################

qiime tools export --input-path demux-single-end.qza --output-path before_cutadapt
qiime tools export --input-path trimmed.qza           --output-path after_cutadapt

conda activate base

# Stats and primer counts as above

cat before_cutadapt/*.fastq.gz > before_cutadapt/combined.fastq.gz
seqkit stats before_cutadapt/combined.fastq.gz >> stats_primers.txt
rm before_cutadapt/combined.fastq.gz

cat after_cutadapt/*.fastq.gz > after_cutadapt/combined.fastq.gz
seqkit stats after_cutadapt/combined.fastq.gz >> stats_primers.txt
rm after_cutadapt/combined.fastq.gz

cat before_cutadapt/*.fastq.gz | \
  seqkit grep -s -i -d -p CTTGGTCATTTAGAGGAAGTAA -C >> count_primers.txt
cat after_cutadapt/*.fastq.gz | \
  seqkit grep -s -i -d -p CTTGGTCATTTAGAGGAAGTAA -C >> count_primers.txt

cat before_cutadapt/*.fastq.gz | \
  seqkit grep -s -i -d -p GCATCGATGAAGAACGCAGC -C >> count_primers.txt
cat after_cutadapt/*.fastq.gz | \
  seqkit grep -s -i -d -p GCATCGATGAAGAACGCAGC -C >> count_primers.txt

# 3. Denoising / ASV inference (DADA2)
# ------------------------------------

qiime dada2 denoise-single \
  --i-demultiplexed-seqs trimmed.qza \
  --p-trunc-len 0 \
  --p-trim-left 0 \
  --p-trunc-q 8 \
  --p-max-ee 8.0 \
  --p-n-threads 18 \
  --p-pooling-method pseudo \
  --output-dir dada2 \
  --verbose

qiime tools export \
  --input-path dada2/denoising_stats.qza \
  --output-path dada2

qiime metadata tabulate \
  --m-input-file dada2/denoising_stats.qza \
  --o-visualization dada2/stats-dada2.qzv

qiime feature-table summarize \
  --i-table dada2/table.qza \
  --o-visualization dada2/table.qzv

cd ..

#################### MERGE RUNS #####################################

# NOTE: We do not run decontam here due to insufficient negative controls.
# DNA-concentration based decontam was not convincing for this project.

mkdir -p dada2

# Merge feature tables
qiime feature-table merge \
  --i-tables run1/dada2/table.qza \
  --i-tables run2/dada2/table.qza \
  --o-merged-table dada2/table.qza

# Merge representative sequences
qiime feature-table merge-seqs \
  --i-data run1/dada2/representative_sequences.qza \
  --i-data run2/dada2/representative_sequences.qza \
  --o-merged-data dada2/representative_sequences.qza

# Summaries
qiime feature-table summarize \
  --i-table dada2/table.qza \
  --o-visualization dada2/table.qzv

qiime feature-table tabulate-seqs \
  --i-data dada2/representative_sequences.qza \
  --o-visualization dada2/representative_sequences.qzv

# qiime tools view dada2/table.qzv
# qiime tools view dada2/representative_sequences.qzv


#################### TAXONOMIC ASSIGNMENT ###########################

# UNITE 9.0 – dynamic, all taxa 

mkdir -p taxo

qiime feature-classifier classify-sklearn \
  --i-classifier ~/qiime2/analysis/databases/ITS/UNITE9.0_q2024.5/unite_ver9_dynamic_all-Q2-2024.5.qza \
  --i-reads dada2/representative_sequences.qza \
  --o-classification taxo/taxonomy.qza \
  --p-n-jobs 1

qiime metadata tabulate \
  --m-input-file taxo/taxonomy.qza \
  --o-visualization taxo/taxonomy.qzv

qiime taxa barplot \
  --i-table dada2/table.qza \
  --i-taxonomy taxo/taxonomy.qza \
  --m-metadata-file ../metadata_all.tsv \
  --o-visualization taxo/taxa_barplot_all.qzv


#################### FILTER ASVs ####################################

mkdir -p filterASV

# Filter out non-fungal ASVs (keep only Fungi)
qiime taxa filter-table \
  --i-table dada2/table.qza \
  --i-taxonomy taxo/taxonomy.qza \
  --o-filtered-table filterASV/taxfiltered_F_table.qza \
  --p-include Fungi

# Further filter to keep ASVs classified at phylum level (prefix "p_")
qiime taxa filter-table \
  --i-table filterASV/taxfiltered_F_table.qza \
  --i-taxonomy taxo/taxonomy.qza \
  --o-filtered-table filterASV/taxfiltered_table.qza \
  --p-include p_

qiime feature-table summarize \
  --i-table filterASV/taxfiltered_table.qza \
  --o-visualization filterASV/taxfiltered_table.qzv

# Filter representative sequences accordingly
qiime feature-table filter-seqs \
  --i-data dada2/representative_sequences.qza \
  --i-table filterASV/taxfiltered_table.qza \
  --o-filtered-data filterASV/taxfiltered_rep_seq.qza

qiime feature-table tabulate-seqs \
  --i-data filterASV/taxfiltered_rep_seq.qza \
  --o-visualization filterASV/taxfiltered_rep_seq.qzv


#################### EXPORT FOR qiime2R / R ANALYSIS ################

mkdir -p qiime2R_ITS

# Copy DADA2 results
cp -R dada2    qiime2R_ITS/dada2

# Copy filtered ASVs
cp -R filterASV qiime2R_ITS/filterASV

# Copy taxonomy
cp -R taxo        qiime2R_ITS/taxo_UNITE10_all
cp -R taxo10nolyF qiime2R_ITS/taxo_UNITE10_Fungi

# Export representative sequences as FASTA
qiime tools export \
  --input-path qiime2R_ITS/dada2/representative_sequences.qza \
  --output-path qiime2R_ITS/dada2/rep_seqs

mv qiime2R_ITS/dada2/rep_seqs/dna-sequences.fasta \
   qiime2R_ITS/dada2/ASV-seqs.fasta

rm -R qiime2R_ITS/dada2/rep_seqs

# Metadata (shared across runs)
cp ../metadata_all.tsv qiime2R_ITS/metadata.tsv

# Compress for transfer to R
zip -r qiime2R_ITS_forR ./qiime2R_ITS
