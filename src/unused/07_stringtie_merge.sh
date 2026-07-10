#!/usr/bin/env bash
set -Eeuo pipefail
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Load workflow configuration.
source "${SCRIPT_DIR}/01_common.sh"
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate "${ENV}"
log "Running ${0##*/}"

gtflist="${COMBINED_DIR}/gtf_list.txt"
n_gtf="$(wc -l < "${gtflist}")"
if (( n_gtf == 0 )); then
    echo "No per-sample GTFs were produced" >&2
    exit 1
fi
log "Per-sample GTFs ready: ${n_gtf}"

stringtie_args=(
    -p "${THREADS}"
    -F 0
    -g 0 
    -f 0 
    -i 
    -m 0 
    -T 0 
    -c 0
)

# --- 2. stringtie --merge with GENCODE as guide ---------------------------
# CRITICAL: BED-derived GTFs have no FPKM/TPM/coverage, so the default merge
# thresholds (-F 1.0, -T 1.0) would discard every transcript. Zero them all.
# -i keeps retained-intron isoforms; -m 0 imposes no length floor.
log "Running stringtie --merge"
stringtie --merge \
    "${stringtie_args[@]}" \
    -G "${GTF}" \
    -o "${MERGED_GTF}" \
    "${gtflist}"

n_tx="$(awk -F'\t' '$3=="transcript"' "${MERGED_GTF}" | wc -l)"
n_ref="$(awk -F'\t' '$3=="transcript"' "${GTF}" | wc -l)"
log "Merged GTF: ${MERGED_GTF}"
log "Transcripts in merge: ${n_tx}  (GENCODE reference has ${n_ref})"
