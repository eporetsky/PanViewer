#!/usr/bin/env bash
# List all directed BLAST tasks including self: QUERY\tSUBJECT (N² lines).
# Usage: genetribe_build_blast_tasks.sh <analysis_dir> > blast_tasks.tsv
set -euo pipefail

ANALYSIS="${1:?usage: $0 <analysis_dir>}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LIST="${REPO_ROOT}/scripts/genetribe_list_accessions.sh"

mapfile -t accs < <(bash "${LIST}" "${ANALYSIS}")
n=${#accs[@]}
[[ "${n}" -ge 1 ]] || {
  echo "ERROR: no accessions found" >&2
  exit 1
}

for q in "${accs[@]}"; do
  for s in "${accs[@]}"; do
    printf '%s\t%s\n' "${q}" "${s}"
  done
done
