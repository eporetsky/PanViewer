#!/usr/bin/env bash
# Build per-accession longest-protein FASTA + BLAST DB for GeneTribe -d reuse.
#
# Layout:
#   <analysis>/work/blast/<acc>_long.fa
#   <analysis>/work/blast/db/<acc>.*   (makeblastdb)
#
# Usage: genetribe_prep_blast_dbs.sh <analysis_dir>
# Env:   WORK=  GENETRIBE_EXTRA='-s @'  (separator must match pair jobs)
set -euo pipefail

ANALYSIS=$(realpath "${1:?usage: $0 <analysis_dir>}")
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK=$(realpath -m "${WORK:-${ANALYSIS}/work}")
BLAST_DIR="${GENETRIBE_BLAST_DIR:-${WORK}/blast}"
DB_DIR="${BLAST_DIR}/db"
LIST="${REPO_ROOT}/scripts/genetribe_list_accessions.sh"
SEP_HELPER="${REPO_ROOT}/scripts/genetribe_blast_sep.sh"

mkdir -p "${BLAST_DIR}" "${DB_DIR}"

# Prefer genetribe-upstream longestfasta (same as core)
LONGFASTA=""
for cand in \
  "${GENETRIBE_HOME:-}/bin/longestfasta" \
  "${REPO_ROOT}/genetribe-upstream/bin/longestfasta" \
  "$(command -v longestfasta 2>/dev/null || true)"
do
  if [[ -n "${cand}" && -x "${cand}" ]]; then
    LONGFASTA="${cand}"
    break
  elif [[ -n "${cand}" && -f "${cand}" ]]; then
    LONGFASTA="${cand}"
    break
  fi
done
[[ -n "${LONGFASTA}" ]] || {
  echo "ERROR: longestfasta not found. Run: make setup" >&2
  exit 1
}
command -v makeblastdb >/dev/null 2>&1 || {
  echo "ERROR: makeblastdb not on PATH (activate genetribe conda env)" >&2
  exit 1
}

SEP=$(bash "${SEP_HELPER}")
echo "Preparing BLAST DBs under ${BLAST_DIR} (separator=${SEP})"

n=0
while IFS= read -r acc; do
  [[ -z "${acc}" ]] && continue
  fa="${ANALYSIS}/accessions/${acc}/${acc}.fa"
  [[ -f "${fa}" ]] || {
    echo "ERROR: missing ${fa}" >&2
    exit 1
  }
  long_fa="${BLAST_DIR}/${acc}_long.fa"
  db_prefix="${DB_DIR}/${acc}"

  if [[ ! -f "${long_fa}" || "${fa}" -nt "${long_fa}" ]]; then
    echo "  longestfasta ${acc}"
    python3 "${LONGFASTA}" -i "${fa}" -s "${SEP}" >"${long_fa}.tmp"
    mv "${long_fa}.tmp" "${long_fa}"
  fi

  # Rebuild DB if missing or stale vs long.fa
  if [[ ! -f "${db_prefix}.psq" || "${long_fa}" -nt "${db_prefix}.psq" ]]; then
    echo "  makeblastdb ${acc}"
    makeblastdb -in "${long_fa}" -parse_seqids -hash_index -dbtype prot -out "${db_prefix}" \
      >/dev/null
  fi
  # GeneTribe also accepts <acc>_long.fa next to inputs; keep a stable copy name
  n=$((n + 1))
done < <(bash "${LIST}" "${ANALYSIS}")

echo "Prepared ${n} BLAST DB(s) in ${BLAST_DIR}"
