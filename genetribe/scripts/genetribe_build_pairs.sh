#!/usr/bin/env bash
# Build unordered accession-pair list (A\tB with A < B lexicographically).
# Usage: genetribe_build_pairs.sh <analysis_dir> > pairs.tsv
set -euo pipefail

ANALYSIS="${1:?usage: $0 <analysis_dir>}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LIST="${REPO_ROOT}/scripts/genetribe_list_accessions.sh"

mapfile -t accs < <(bash "${LIST}" "${ANALYSIS}")
n=${#accs[@]}
[[ "${n}" -ge 2 ]] || {
  echo "ERROR: need at least 2 accessions (found ${n})" >&2
  exit 1
}

for ((i = 0; i < n; i++)); do
  for ((j = i + 1; j < n; j++)); do
    a="${accs[i]}"
    b="${accs[j]}"
    if [[ "${a}" < "${b}" ]]; then
      printf '%s\t%s\n' "${a}" "${b}"
    else
      printf '%s\t%s\n' "${b}" "${a}"
    fi
  done
done
