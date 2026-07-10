#!/usr/bin/env bash
set -Eeuo pipefail
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Load workflow configuration.
source "${SCRIPT_DIR}/01_common.sh"
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate "${ENV}"
log "Running ${0##*/}"

ps_dir="${COMBINED_DIR}/persample_gtf"
mkdir -p "${ps_dir}"

# --- 1. Convert each per-sample lncRNA BED -> GTF (exon structure kept) ----
# Per-sample MSTRG ids collide across files; that's fine, stringtie --merge
# reconciles by structure, not by id (it's the normal multi-sample case).
shopt -s nullglob
beds=("${LNCRNA_BED_DIR}"/*_true_lncRNA.rf.bed)
shopt -u nullglob
if (( ${#beds[@]} == 0 )); then
    echo "No lncRNA BED files found under ${LNCRNA_BED_DIR}" >&2
    exit 1
fi
log "Converting ${#beds[@]} per-sample BEDs to GTF"

gtflist="${COMBINED_DIR}/gtf_list.txt"
: > "${gtflist}"

for bed in "${beds[@]}"; do
    s="$(basename "${bed}" .bed)"
    tmp_bed="${ps_dir}/${s}.bed"
    gp="${ps_dir}/${s}.genePred"
    gtf="${ps_dir}/${s}.gtf"

    # Strip header/track lines; null thick (cols 7,8) so transcripts convert
    # as noncoding; keep BED12 block columns intact for exon structure.
    awk -F'\t' 'BEGIN{OFS="\t"}
        $1!="chrom" && $4!="name" && $1 !~ /^(track|browser|#)/ && NF>=12 {
            $7=$2; $8=$2; print
        }' "${bed}" > "${tmp_bed}"

    if [[ ! -s "${tmp_bed}" ]]; then
        log "WARNING: no usable BED12 rows in ${s}, skipping"
        rm -f "${tmp_bed}"
        continue
    fi

    bedToGenePred "${tmp_bed}" "${gp}"
    genePredToGtf -source=Flnc file "${gp}" "${gtf}"
    rm -f "${tmp_bed}" "${gp}"
    printf '%s\n' "${gtf}" >> "${gtflist}"
done

n_gtf="$(wc -l < "${gtflist}")"
if (( n_gtf == 0 )); then
    echo "No per-sample GTFs were produced" >&2
    exit 1
fi
log "Per-sample GTFs ready: ${n_gtf}"
