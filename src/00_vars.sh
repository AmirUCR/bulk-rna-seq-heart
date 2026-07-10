#!/usr/bin/env bash

# Central configuration for the RNA-seq workflow.
# This file is meant to be sourced by the other workflow scripts:
#
#   source "$(dirname "$0")/00_vars.sh"
#
# Note:
# - Paths are resolved relative to this file, not the caller's current directory.

export WORKFLOW_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PROJECT_DIR="$(cd "${WORKFLOW_DIR}/.." && pwd)"

export ENV="heart_env"
export THREADS=64
export STRAND_THRESH=0.75

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
export DATASET="hg38"
export REF_EXTENSION=".fa"
export REF_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
export GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
export CHROM_SIZES_URL="https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes"
export REF="${GENOMIC}/${DATASET}${REF_EXTENSION}"
export GTF="${GENOMIC}/gencode.v49.annotation.gtf"
export CHROM_SIZES="${GENOMIC}/${DATASET}.chrom.sizes"

# OUT
export IDX_DIR="${GENOMIC}/hisat2_index"
export HISAT2_IDX="${IDX_DIR}/${DATASET}_index"
export COMBINED_DIR="${LOCAL_RESULTS_DIR}/combined"

export COUNTS_OUT="${LOCAL_RESULTS_DIR}/counts"
export WILCOXON_OUT="${LOCAL_RESULTS_DIR}/wilcoxon_de"
