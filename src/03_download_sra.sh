#!/usr/bin/env bash
set -Eeuo pipefail
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Load workflow configuration first.
source "${SCRIPT_DIR}/01_common.sh"
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate "${ENV}"
log "Running ${0##*/}"

# Create an output directory for FASTQs
mkdir -p "${READS_UNTRIM}"

while IFS= read -r ACC || [[ -n "${ACC}" ]]; do
    # Skip empty lines or comments
    [[ -z "${ACC}" || "${ACC}" =~ ^# ]] && continue

    log "------------------------------------------------------"
    log "Processing: ${ACC}"
    log "------------------------------------------------------"

    # 1. Prefetch the data
    log "Starting prefetch..."
    prefetch "${ACC}" --progress

    # 2. Extract to FASTQ
    # Using --split-3 to handle paired-end and single-end automatically
    if fasterq-dump --split-3 --threads ${THREADS} --outdir "${READS_UNTRIM}" "${ACC}"; then
        log "Compressing files..."
        pigz "${READS_UNTRIM}"/"${ACC}"*.fastq
        rm -f "${ACC}/${ACC}.sra"
        rmdir "${ACC}" 2>/dev/null || true
    else
        log "FAILED: ${ACC}" >> "failed_${ACC}.txt"
    fi

    log "Done with ${ACC}"
done < "${ACCESSIONS_FILE}"

log "------------------------------------------------------"
log "Batch processing complete. Files are in ${READS_UNTRIM}"
