#!/usr/bin/env bash

# Central configuration for the RNA-seq workflow.
# This file is meant to be sourced by the other workflow scripts:
#
#   source "$(dirname "$0")/00_vars.sh"
#
# Note:
# - Paths are resolved relative to this file, not the caller's current directory.

export THREADS=64
export ENV="heart_env"
export REF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/GRCh38.primary_assembly.genome.fa.gz"
export GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.primary_assembly.annotation.gtf.gz"
# export CHROM_SIZES_URL="https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes"

# Used in 09_find_strandedness
export STRAND_THRESH=0.75

# Used in 06_trim.sh
export FASTP_ARGS=(
    --length_required 25
    --trim_poly_x
    --poly_x_min_len 10
)

# Used in 12_featurecounts.sh
export FEATURECOUNTS_ARGS=(
    -Q 10
)

# --- PATHS ---
export WORKFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJECT_DIR="$(cd "${WORKFLOW_DIR}/.." && pwd)"

export LOCAL_DATA_DIR="${WORKFLOW_DIR}/data"
export SHARED_DATA_DIR="${PROJECT_DIR}/data"

export READS_TRIM="${SHARED_DATA_DIR}/reads/trimmed"
export READS_UNTRIM="${SHARED_DATA_DIR}/reads/untrimmed"

export SHARED_BAM_DIR="${PROJECT_DIR}/results/bam"

export LOCAL_RESULTS_DIR="${WORKFLOW_DIR}/results"
export READS_OUT_UNTRIM="${LOCAL_RESULTS_DIR}/reads/untrimmed"
export READS_OUT_TRIM="${LOCAL_RESULTS_DIR}/reads/trimmed"
export STRAND_DIR="${LOCAL_RESULTS_DIR}/trimmed/strand"

export ACCESSIONS_FILE="${SHARED_DATA_DIR}/accessions.txt"
export SAMPLES_FILE="${SHARED_DATA_DIR}/samples.tsv"

# REF
export GENOMIC="${SHARED_DATA_DIR}/genomic"
export DATASET="$(basename "${REF_URL}" .fa.gz)"
export REF="${GENOMIC}/$(basename "${REF_URL}" .gz)"
export GTF="${GENOMIC}/$(basename "${GTF_URL}" .gz)"
# export CHROM_SIZES="${GENOMIC}/$(basename "${CHROM_SIZES_URL}")"

# OUT
export IDX_DIR="${GENOMIC}/hisat2_index"
export HISAT2_IDX="${IDX_DIR}/${DATASET}_index"
export COMBINED_DIR="${LOCAL_RESULTS_DIR}/combined"

export COUNTS_OUT="${LOCAL_RESULTS_DIR}/counts"
export WILCOXON_OUT="${LOCAL_RESULTS_DIR}/wilcoxon_de"
