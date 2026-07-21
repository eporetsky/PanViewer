#!/usr/bin/env bash
# Rewrite existing accessions/<acc>/<acc>.bed to GeneTribe 6-col format in place.
# Use after a prep that left 7-col Pandagma BED symlinks/files.
#
# Usage: genetribe_fix_beds.sh <analysis_dir>
set -euo pipefail

ANALYSIS=$(realpath "${1:?usage: $0 <analysis_dir>}")
ACC_ROOT="${ANALYSIS}/accessions"
[[ -d "${ACC_ROOT}" ]] || { echo "ERROR: missing ${ACC_ROOT}" >&2; exit 1; }

n=0
fixed=0
while IFS= read -r -d '' bed; do
  n=$((n + 1))
  acc=$(basename "$(dirname "${bed}")")
  ncols=$(gawk 'NF && $1 !~ /^#/ { print NF; exit }' "${bed}" || true)
  if [[ "${ncols}" == "6" ]]; then
    continue
  fi
  if [[ "${ncols}" != "7" ]]; then
    echo "ERROR: ${bed} has ${ncols:-0} columns (need 6 or 7)" >&2
    exit 1
  fi
  tmp="${bed}.gt6.$$"
  # Resolve symlink source if needed, then rewrite as a real file
  src="${bed}"
  if [[ -L "${bed}" ]]; then
    src=$(readlink -f "${bed}")
  fi
  gawk -vOFS='\t' '
    /^#/ || NF == 0 { next }
    NF >= 7 { print $1, $2, $3, $7, $5, $6; next }
    { print }
  ' "${src}" >"${tmp}"
  rm -f "${bed}"
  mv "${tmp}" "${bed}"
  echo "fixed ${acc}: 7 → 6 columns"
  fixed=$((fixed + 1))
done < <(find "${ACC_ROOT}" -mindepth 2 -maxdepth 2 -type f -name '*.bed' -print0 -o -type l -name '*.bed' -print0)

echo "Checked ${n} BED(s); converted ${fixed} to 6-col GeneTribe format."
