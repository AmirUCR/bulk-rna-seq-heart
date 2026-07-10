#!/usr/bin/env bash
set -Eeuo pipefail
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Load workflow configuration.
source "${SCRIPT_DIR}/01_common.sh"
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate "${ENV}"
log "Running ${0##*/}"

# Discover the trimmed fastq.gz files under READS_TRIM and group them into
# samples. A *_1 / *_2 pair is paired-end; everything else is single-end.
shopt -s nullglob
fastq_files=("${READS_TRIM}"/*.fastq.gz)
shopt -u nullglob

declare -a JOBS=()
for f in "${fastq_files[@]}"; do
    base="$(basename "${f}")"
    case "${base}" in
        *_2.fastq.gz)
            # _2 of a pair: handled together with its _1 mate, so skip here.
            continue
            ;;
        *_1.fastq.gz)
            sample_id="${base%_1.fastq.gz}"
            r2="${READS_TRIM}/${sample_id}_2.fastq.gz"
            if [[ -f "${r2}" ]]; then
                JOBS+=( "paired"$'\t'"${sample_id}"$'\t'"${f}"$'\t'"${r2}" )
            else
                log "WARNING: ${base} has no _2 mate, treating as single-end"
                JOBS+=( "single"$'\t'"${sample_id}"$'\t'"${f}" )
            fi
            ;;
        *.fastq.gz)
            sample_id="${base%.fastq.gz}"
            JOBS+=( "single"$'\t'"${sample_id}"$'\t'"${f}" )
            ;;
    esac
done

NUM_SAMPLES="${#JOBS[@]}"
if (( NUM_SAMPLES == 0 )); then
    log "No samples found in ${READS_TRIM}" >&2
    exit 1
fi

# Create per-sample output directories.
for job in "${JOBS[@]}"; do
    IFS=$'\t' read -r layout sample_id r1 r2 <<< "${job}"
    mkdir -p "${SHARED_BAM_DIR}/${sample_id}"
done

infer_hisat2_strand
log "Detected hisat2 strandedness: ${HISAT2_STRANDEDNESS}"

align() {
    local layout="$1"
    local r1="$2"
    local r2="$3"   # empty for single-end
    local tag="$4"

    # hisat2 strandedness: two-letter for paired (FR/RF), single-letter for
    # unpaired (F/R), and omitted entirely when unstranded ("none").
    local strand_args=()
    if [[ "${HISAT2_STRANDEDNESS}" != "none" ]]; then
        if [[ "${layout}" == "paired" ]]; then
            strand_args=(--rna-strandness "${HISAT2_STRANDEDNESS}")
        else
            strand_args=(--rna-strandness "${HISAT2_STRANDEDNESS:0:1}")
        fi
    fi

    if [[ "${layout}" == "paired" ]]; then
        hisat2 \
            -p "${THREADS}" \
            "${strand_args[@]}" \
            -x "${HISAT2_IDX}" \
            -1 "${r1}" \
            -2 "${r2}" \
            2> "${tag}.hisat2.log" \
            | samtools sort -@ "${THREADS}" -o "${tag}.sorted.bam" -
    else
        hisat2 \
            -p "${THREADS}" \
            "${strand_args[@]}" \
            -x "${HISAT2_IDX}" \
            -U "${r1}" \
            2> "${tag}.hisat2.log" \
            | samtools sort -@ "${THREADS}" -o "${tag}.sorted.bam" -
    fi

    samtools index "${tag}.sorted.bam"
}

# —— Run
for job in "${JOBS[@]}"; do
    IFS=$'\t' read -r layout sample_id r1 r2 <<< "${job}"
    align \
        "${layout}" \
        "${r1}" \
        "${r2}" \
        "${SHARED_BAM_DIR}/${sample_id}/trimmed"
done
