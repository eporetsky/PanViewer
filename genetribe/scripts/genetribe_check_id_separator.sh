#!/usr/bin/env bash
# Require gene-level IDs: BED col4 must match FASTA headers exactly.
# This workflow does NOT use GeneTribe transcript stripping (we always pass
# -s @ so "." inside names like TraesPARA_EIv1.0_* is left alone).
#
# Usage:
#   genetribe_check_id_separator.sh <analysis_dir>
#   genetribe_check_id_separator.sh <acc.fa> <acc.bed>
#
# Env:
#   GENETRIBE_SKIP_ID_CHECK=1   skip (not recommended)
set -euo pipefail

if [[ "${GENETRIBE_SKIP_ID_CHECK:-0}" == "1" ]]; then
  echo "WARN: GENETRIBE_SKIP_ID_CHECK=1 — skipping gene-ID check" >&2
  exit 0
fi

check_one() {
  local fa="$1" bed="$2" label="${3:-$1}"
  local tmp
  tmp=$(mktemp -d)
  awk -F'\t' 'NF>=4 && $1 !~ /^#/ { print $4 }' "${bed}" | sort -u >"${tmp}/bed_ids"
  grep '^>' "${fa}" | sed 's/^>//;s/[[:space:]].*//' | sort -u >"${tmp}/fa_ids"

  local n_bed n_fa n_both
  n_bed=$(wc -l <"${tmp}/bed_ids" | tr -d ' ')
  n_fa=$(wc -l <"${tmp}/fa_ids" | tr -d ' ')
  n_both=$(comm -12 "${tmp}/bed_ids" "${tmp}/fa_ids" | wc -l | tr -d ' ')
  rm -rf "${tmp}"

  if [[ "${n_bed}" -eq 0 || "${n_fa}" -eq 0 ]]; then
    echo "ERROR: ${label}: empty BED gene IDs (${n_bed}) or FASTA IDs (${n_fa})" >&2
    return 1
  fi
  if [[ "${n_both}" -eq 0 ]]; then
    echo "ERROR: ${label}: BED col4 and FASTA headers share 0 IDs." >&2
    echo "  Provide gene-level IDs only (same string in .bed col4 and FASTA header)." >&2
    echo "  Do not rely on GeneTribe -s transcript stripping in this workflow." >&2
    echo "  bed=${bed}" >&2
    echo "  fa=${fa}" >&2
    return 1
  fi
  # Soft warning if overlap is tiny relative to either set
  if [[ "${n_both}" -lt $((n_bed / 2)) || "${n_both}" -lt $((n_fa / 2)) ]]; then
    echo "WARN: ${label}: only ${n_both} shared IDs (bed=${n_bed} fa=${n_fa}); check isoform stripping." >&2
  fi
  echo "OK ${label}: shared_ids=${n_both} bed=${n_bed} fa=${n_fa}"
  return 0
}

if [[ $# -eq 1 ]]; then
  ANALYSIS=$(realpath "$1")
  ACC="${ANALYSIS}/accessions"
  [[ -d "${ACC}" ]] || { echo "ERROR: missing ${ACC}" >&2; exit 1; }
  echo "Checking gene-level BED↔FASTA ID match under ${ACC}"
  fail=0
  while IFS= read -r -d '' bed; do
    acc=$(basename "$(dirname "${bed}")")
    fa="${ACC}/${acc}/${acc}.fa"
    [[ -f "${fa}" ]] || { echo "ERROR: missing ${fa}" >&2; fail=1; continue; }
    check_one "${fa}" "${bed}" "${acc}" || fail=1
  done < <(find "${ACC}" -mindepth 2 -maxdepth 2 \( -type f -o -type l \) -name '*.bed' -print0 | sort -z)
  [[ "${fail}" -eq 0 ]] || exit 1
  echo "All accessions OK (gene-level IDs)."
  exit 0
elif [[ $# -eq 2 ]]; then
  check_one "$(realpath "$1")" "$(realpath "$2")" "$(basename "$2" .bed)"
else
  echo "usage: $0 <analysis_dir> | $0 <acc.fa> <acc.bed>" >&2
  exit 1
fi
