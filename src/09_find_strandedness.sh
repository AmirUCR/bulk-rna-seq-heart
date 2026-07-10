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

if (( ${#JOBS[@]} == 0 )); then
    echo "No samples found in ${READS_TRIM}" >&2
    exit 1
fi

process_sample() {
    local layout="$1"
    local sample_id="$2"
    local r1="$3"
    local r2="${4:-}"   # empty for single-end
    local sdir="${STRAND_DIR}/${sample_id}"

    # Subsample reads (reproducible, random selection) and align.
    if [[ "${layout}" == "paired" ]]; then
        log "Sampling reads for ${sample_id} (paired-end)"
        seqtk sample -s100 "${r1}" 500000 > "${sdir}/1_500k_sample.fq"
        seqtk sample -s100 "${r2}" 500000 > "${sdir}/2_500k_sample.fq"
        log "Aligning ${sample_id}"
        hisat2 -p "${THREADS}" -x "${HISAT2_IDX}" \
            -1 "${sdir}/1_500k_sample.fq" \
            -2 "${sdir}/2_500k_sample.fq" \
            -S "${sdir}/sample.sam"
        rm -f "${sdir}/1_500k_sample.fq" "${sdir}/2_500k_sample.fq"
    else
        log "Sampling reads for ${sample_id} (single-end)"
        seqtk sample -s100 "${r1}" 500000 > "${sdir}/reads_500k_sample.fq"
        log "Aligning ${sample_id}"
        hisat2 -p "${THREADS}" -x "${HISAT2_IDX}" \
            -U "${sdir}/reads_500k_sample.fq" \
            -S "${sdir}/sample.sam"
        rm -f "${sdir}/reads_500k_sample.fq"
    fi

    samtools view -bS "${sdir}/sample.sam" | \
        samtools sort -@ "${THREADS}" -m 2G -T "${sdir}/tmp.sort" \
            -o "${sdir}/sample.sorted.bam"
    rm -f "${sdir}/sample.sam"
    samtools index "${sdir}/sample.sorted.bam"

    # Run infer_experiment.py
    log "Running infer_experiment for ${sample_id}"
    infer_experiment.py \
        -i "${sdir}/sample.sorted.bam" \
        -r "${GENOMIC}/${DATASET}.bed" > "${sdir}/${sample_id}_strandedness.txt"
    rm -f "${sdir}/sample.sorted.bam" "${sdir}/sample.sorted.bam.bai"
}

# Prepare gene annotation in BED12 format (UCSC tools, required by RSeQC).
log "Creating ${GENOMIC}/${DATASET}.bed"
awk '$3 != "gene"' "${GTF}" | grep -v '^#' \
    | gtfToGenePred /dev/stdin /dev/stdout \
    | genePredToBed /dev/stdin "${GENOMIC}/${DATASET}.bed"

# Strandedness is a library-prep property, so infer it from a single
# representative sample rather than all of them.
IFS=$'\t' read -r layout sample_id r1 r2 <<< "${JOBS[0]}"
log "Inferring strandedness from one sample: ${sample_id}"
mkdir -p "${STRAND_DIR}/${sample_id}"
process_sample "${layout}" "${sample_id}" "${r1}" "${r2}"

log "Done. Result: ${STRAND_DIR}/${sample_id}/${sample_id}_strandedness.txt"
