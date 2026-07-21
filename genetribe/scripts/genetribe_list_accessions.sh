#!/usr/bin/env bash
# List accession prefixes under <analysis>/accessions/<acc>/<acc>.{fa,bed,chrlist}.
# Prints one accession id (basename) per line, sorted.
#
# Usage: genetribe_list_accessions.sh <analysis_dir>
set -euo pipefail

ANALYSIS="${1:?usage: $0 <analysis_dir>}"
ACC_DIR="${ANALYSIS}/accessions"

[[ -d "${ACC_DIR}" ]] || {
  echo "ERROR: missing ${ACC_DIR}" >&2
  exit 1
}

found=0
while IFS= read -r -d '' d; do
  acc=$(basename "${d}")
  fa="${d}/${acc}.fa"
  bed="${d}/${acc}.bed"
  chrlist="${d}/${acc}.chrlist"
  if [[ ! -f "${fa}" || ! -f "${bed}" || ! -f "${chrlist}" ]]; then
    echo "ERROR: accession '${acc}' needs ${acc}.fa, ${acc}.bed, and ${acc}.chrlist under ${d}" >&2
    exit 1
  fi
  printf '%s\n' "${acc}"
  found=1
done < <(find "${ACC_DIR}" -mindepth 1 -maxdepth 1 -type d -print0 | sort -z)

[[ "${found}" -eq 1 ]] || {
  echo "ERROR: no accession directories under ${ACC_DIR}" >&2
  exit 1
}
