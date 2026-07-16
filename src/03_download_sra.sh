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

while IFS= read -r acc || [[ -n "${acc}" ]]; do
    # Skip empty lines or comments
    [[ -z "${acc}" || "${acc}" =~ ^# ]] && continue

    log "------------------------------------------------------"
    log "Processing: ${acc}"
    log "------------------------------------------------------"

    # 1. Prefetch the data
    log "Starting prefetch..."
    prefetch "${acc}" --progress

    # 2. Extract to FASTQ
    # Using --split-3 to handle paired-end and single-end automatically
    if fasterq-dump --split-3 --threads ${THREADS} --outdir "${READS_UNTRIM}" "${acc}"; then
        log "Compressing files..."
        pigz "${READS_UNTRIM}"/"${acc}"*.fastq
        rm -f "${acc}/${acc}.sra"
        rmdir "${acc}" 2>/dev/null || true
    else
        log "FAILED: ${acc}" >> "failed_${acc}.txt"
    fi

    log "Done with ${acc}"
done < "${ACCESSIONS_FILE}"

log "------------------------------------------------------"
log "Batch processing complete. Files are in ${READS_UNTRIM}"

log "Done."
