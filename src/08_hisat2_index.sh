#!/usr/bin/env bash
set -Eeuo pipefail
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Load workflow configuration first.
source "${SCRIPT_DIR}/01_common.sh"
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate "${ENV}"
log "Running ${0##*/}"

mkdir -p "${IDX_DIR}"

# —— Build HISAT2 index (only if missing)
if ! ls "${HISAT2_IDX}".* >/dev/null 2>&1; then
  log "Running hisat2-build"
  
  hisat2-build -p "${THREADS}" "${REF}" "${HISAT2_IDX}"
fi
