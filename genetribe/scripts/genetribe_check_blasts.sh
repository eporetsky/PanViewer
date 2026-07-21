#!/usr/bin/env bash
# Verify all directed BLAST files exist for GeneTribe -d (N² including self).
# Usage: genetribe_check_blasts.sh <analysis_dir>
# Exit 0 if complete; exit 1 and list missing count otherwise.
set -euo pipefail

ANALYSIS=$(realpath "${1:?usage: $0 <analysis_dir>}")
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK=$(realpath -m "${WORK:-${ANALYSIS}/work}")
BLAST_DIR=$(realpath -m "${GENETRIBE_BLAST_DIR:-${WORK}/blast}")
TASKS="${REPO_ROOT}/scripts/genetribe_build_blast_tasks.sh"

[[ -d "${BLAST_DIR}" ]] || {
  echo "ERROR: missing blast dir ${BLAST_DIR}" >&2
  exit 1
}

missing=0
total=0
while IFS=$'\t' read -r q s; do
  [[ -z "${q:-}" || -z "${s:-}" ]] && continue
  total=$((total + 1))
  f="${BLAST_DIR}/${q}_${s}.blast"
  if [[ ! -f "${f}" || ! -s "${f}" ]]; then
    missing=$((missing + 1))
    if [[ "${missing}" -le 20 ]]; then
      echo "MISSING ${q}_${s}.blast" >&2
    fi
  fi
done < <(bash "${TASKS}" "${ANALYSIS}")

have=$((total - missing))
echo "BLAST store ${BLAST_DIR}: ${have}/${total} present"
if [[ "${missing}" -gt 0 ]]; then
  echo "ERROR: ${missing} BLAST file(s) missing (showing up to 20 above)" >&2
  exit 1
fi
echo "BLAST store complete."
