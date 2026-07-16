#!/usr/bin/env bash
set -Eeuo pipefail
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Load workflow configuration first.
source "${SCRIPT_DIR}/01_common.sh"
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate "${ENV}"
log "Running ${0##*/}"

# Create output directory.
mkdir -p "${READS_TRIM}"

# Discover the fastq.gz files under READS_UNTRIM and group them into samples.
# A *_1 / *_2 pair is paired-end; everything else is single-end.
shopt -s nullglob
fastq_files=()
for f in "${READS_UNTRIM}"/*.fastq.gz; do
    base="$(basename "${f}")"
    if [[ -s "${READS_TRIM}/${base}" ]]; then
        log "Skipping ${base} (already trimmed)"
        continue
    fi
    fastq_files+=("${f}")
done
shopt -u nullglob

declare -a jobs=()
for f in "${fastq_files[@]}"; do
    base="$(basename "${f}")"
    case "${base}" in
        *_2.fastq.gz)
            # R2 of a pair: handled together with its _1 mate, so skip here.
            continue
            ;;
        *_1.fastq.gz)
            sample_id="${base%_1.fastq.gz}"
            r2="${READS_UNTRIM}/${sample_id}_2.fastq.gz"
            if [[ -f "${r2}" ]]; then
                jobs+=( "paired"$'\t'"${sample_id}"$'\t'"${f}"$'\t'"${r2}" )
            else
                log "WARNING: ${base} has no _2 mate, treating as single-end"
                jobs+=( "single"$'\t'"${sample_id}"$'\t'"${f}" )
            fi
            ;;
        *.fastq.gz)
            sample_id="${base%.fastq.gz}"
            jobs+=( "single"$'\t'"${sample_id}"$'\t'"${f}" )
            ;;
    esac
done

num_samples="${#jobs[@]}"
if (( num_samples == 0 )); then
    log "No samples found in ${READS_UNTRIM}" >&2
    exit 1
fi

max_jobs="${num_samples}"
if (( max_jobs > THREADS )); then
    max_jobs="${THREADS}"
fi
threads_per_sample=$(( THREADS / max_jobs ))
if (( threads_per_sample < 1 )); then
    threads_per_sample=1
fi

log "Total threads: ${THREADS}"
log "Samples: ${num_samples}"
log "Concurrent fastp jobs: ${max_jobs}"
log "Threads per fastp job: ${threads_per_sample}"

run_fastp_paired() {
    local sample_id="$1"
    local r1="$2"
    local r2="$3"
    log "Running fastp for ${sample_id} (paired-end)"
    fastp \
        -i "${r1}" \
        -I "${r2}" \
        -o "${READS_TRIM}/${sample_id}_1.fastq.gz" \
        -O "${READS_TRIM}/${sample_id}_2.fastq.gz" \
        --html "${READS_TRIM}/fastp_report_${sample_id}.html" \
        --json "${READS_TRIM}/fastp_report_${sample_id}.json" \
        --thread "${threads_per_sample}" \
        "${FASTP_ARGS[@]}"
}

run_fastp_single() {
    local sample_id="$1"
    local r1="$2"
    log "Running fastp for ${sample_id} (single-end)"
    fastp \
        -i "${r1}" \
        -o "${READS_TRIM}/${sample_id}.fastq.gz" \
        --html "${READS_TRIM}/fastp_report_${sample_id}.html" \
        --json "${READS_TRIM}/fastp_report_${sample_id}.json" \
        --thread "${threads_per_sample}" \
        "${FASTP_ARGS[@]}"
}

active_jobs=0
for job in "${jobs[@]}"; do
    IFS=$'\t' read -r layout sample_id r1 r2 <<< "${job}"

    if [[ "${layout}" == "paired" ]]; then
        run_fastp_paired "${sample_id}" "${r1}" "${r2}" &
    else
        run_fastp_single "${sample_id}" "${r1}" &
    fi

    ((active_jobs+=1))
    if (( active_jobs >= max_jobs )); then
        wait
        active_jobs=0
    fi
done
wait

log "Done."
