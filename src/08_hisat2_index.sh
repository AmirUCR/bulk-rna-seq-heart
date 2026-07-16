#!/usr/bin/env bash
set -Eeuo pipefail
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Load workflow configuration first.
source "${SCRIPT_DIR}/01_common.sh"
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate "${ENV}"
log "Running ${0##*/}"

mkdir -p "${IDX_DIR}"

# Extract splice site coordinates if missing
if [ ! -f "${GENOMIC}/genome.ss" ]; then
    log "${GENOMIC}/genome.ss not found. Generating it..."
    hisat2_extract_splice_sites.py "${GTF}" > "${GENOMIC}/genome.ss"
fi

# Extract exon coordinates if missing
if [ ! -f "${GENOMIC}/genome.exon" ]; then
    log "${GENOMIC}/genome.exon not found. Generating it..."
    hisat2_extract_exons.py "${GTF}" > "${GENOMIC}/genome.exon"
fi

# —— Build HISAT2 index only if missing
if ! ls "${HISAT2_IDX}".* >/dev/null 2>&1; then
    log "Running hisat2-build"

    hisat2-build \
        --ss "${GENOMIC}/genome.ss" \
        --exon "${GENOMIC}/genome.exon" \
        -p "${THREADS}" \
        "${REF}" \
        "${HISAT2_IDX}"
fi

log "Done."
