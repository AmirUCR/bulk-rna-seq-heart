#!/usr/bin/env bash
set +o history  # Turn bash history recording off
old_opts=$(set +o)
set -Eeuo pipefail

common_sh_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${common_sh_dir}/00_vars.sh"

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

iterate_samples() {
    tail -n +1 "${SAMPLES_TSV}"
}

require_file() {
    local path="$1"
    if [[ ! -f "${path}" ]]; then
        log "Missing file: ${path}" >&2
        exit 1
    fi
}

require_dir() {
    local path="$1"
    if [[ ! -d "${path}" ]]; then
        log "Missing directory: ${path}" >&2
        exit 1
    fi
}

get_first_strandedness_file() {
    local file=""
    shopt -s nullglob
    for file in "${STRAND_DIR}"/*/*_strandedness.txt; do
        printf '%s\n' "${file}"
        return 0
    done
    log "No strandedness report found under ${STRAND_DIR}" >&2
    return 1
}

parse_strandedness_fractions() {
    local file="$1"
    local forward=""
    local reverse=""
    require_file "${file}"
    # Paired-end infer_experiment labels reads 1++,1--,2+-,2-+ / 1+-,1-+,2++,2--
    # Single-end labels them ++,-- / +-,-+ . Match either form.
    forward="$(grep -E '(1\+\+,1--,2\+-,2-\+|\+\+,--)' "${file}" | awk '{print $NF}')"
    reverse="$(grep -E '(1\+-,1-\+,2\+\+,2--|\+-,-\+)' "${file}" | awk '{print $NF}')"
    if [[ -z "${forward}" || -z "${reverse}" ]]; then
        log "Could not parse strandedness fractions from ${file}" >&2
        return 1
    fi
    printf '%s\t%s\n' "${forward}" "${reverse}"
}

infer_featurecounts_strand() {
    local file=""
    local forward=""
    local reverse=""

    file="$(get_first_strandedness_file)"
    read -r forward reverse < <(parse_strandedness_fractions "${file}")

    if awk -v f="${forward}" -v r="${reverse}" -v t="${STRAND_THRESH}" 'BEGIN { exit !(r > t && r > f) }'; then
        FEATURE_COUNTS_STRANDEDNESS=2
    elif awk -v f="${forward}" -v r="${reverse}" -v t="${STRAND_THRESH}" 'BEGIN { exit !(f > t && f > r) }'; then
        FEATURE_COUNTS_STRANDEDNESS=1
    else
        FEATURE_COUNTS_STRANDEDNESS=0
    fi

    export FEATURE_COUNTS_STRANDEDNESS
}

infer_hisat2_strand() {
    local file=""
    local forward=""
    local reverse=""

    file="$(get_first_strandedness_file)"
    read -r forward reverse < <(parse_strandedness_fractions "${file}")

    if awk -v f="${forward}" -v r="${reverse}" -v t="${STRAND_THRESH}" 'BEGIN { exit !(r > t && r > f) }'; then
        HISAT2_STRANDEDNESS="RF"
    elif awk -v f="${forward}" -v r="${reverse}" -v t="${STRAND_THRESH}" 'BEGIN { exit !(f > t && f > r) }'; then
        HISAT2_STRANDEDNESS="FR"
    else
        HISAT2_STRANDEDNESS="none"
    fi

    export HISAT2_STRANDEDNESS
}
eval "$old_opts"
set -o history  # Turn bash history recording back on
