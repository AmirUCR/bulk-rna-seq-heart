#!/usr/bin/env bash
set -Eeuo pipefail
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Load workflow configuration.
source "${SCRIPT_DIR}/01_common.sh"
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate "${ENV}"
log "Running ${0##*/}"

mkdir -p "${COUNTS_OUT}"

if [[ -s "${COUNTS_OUT}/counts_pe.txt" || -s "${COUNTS_OUT}/counts_se.txt" ]]; then
    log "counts already exist, skipping"
    exit 0
fi

infer_featurecounts_strand
log "Detected featureCounts strandedness: ${FEATURE_COUNTS_STRANDEDNESS}"

log "Collecting bam fles"
shopt -s nullglob
pe_bams=()
se_bams=()
for f in "${SHARED_BAM_DIR}"/*/trimmed.sorted.bam; do
    # >0 reads with the paired flag => paired-end library
    if [[ "$(samtools view -c -f 1 "${f}")" -gt 0 ]]; then
        pe_bams+=("${f}")
    else
        se_bams+=("${f}")
    fi
done
shopt -u nullglob

if (( ${#pe_bams[@]} == 0 && ${#se_bams[@]} == 0 )); then
    echo "No BAMs found under ${SHARED_BAM_DIR}" >&2
    exit 1
fi

common_args=(
    -T "${THREADS}"
    -s "${FEATURE_COUNTS_STRANDEDNESS}"
    -a "${GTF}"
    --extraAttributes "gene_name"
)

if (( ${#pe_bams[@]} > 0 )); then
    log "featureCounts on ${#pe_bams[@]} paired-end BAMs"
    featureCounts \
        "${common_args[@]}" \
        "${FEATURECOUNTS_ARGS[@]}" \
        -p --countReadPairs -C \
        -o "${COUNTS_OUT}/counts_pe.txt" \
        "${pe_bams[@]}"
fi

if (( ${#se_bams[@]} > 0 )); then
    log "featureCounts on ${#se_bams[@]} single-end BAMs"
    featureCounts \
        "${common_args[@]}" \
        "${FEATURECOUNTS_ARGS[@]}" \
        -o "${COUNTS_OUT}/counts_se.txt" \
        "${se_bams[@]}"
fi

log "Done."
