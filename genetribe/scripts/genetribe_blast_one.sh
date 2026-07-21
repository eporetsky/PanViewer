#!/usr/bin/env bash
# Run one directed BLASTP for GeneTribe -d reuse: QUERY vs SUBJECT → QUERY_SUBJECT.blast
#
# Usage: genetribe_blast_one.sh -d <analysis_dir> -q <query_acc> -s <subject_acc>
# Env:   WORK=  GENETRIBE_BLAST_DIR=  GENETRIBE_EVALUE=1e-5  GENETRIBE_THREADS=
#        GENETRIBE_SKIP_EXISTING=1
set -euo pipefail

ANALYSIS=""
Q=""
S=""

usage() {
  echo "usage: $0 -d <analysis_dir> -q <query_acc> -s <subject_acc>" >&2
  exit 1
}

while getopts ":d:q:s:h" opt; do
  case "${opt}" in
    d) ANALYSIS=$(realpath "${OPTARG}") ;;
    q) Q="${OPTARG}" ;;
    s) S="${OPTARG}" ;;
    h) usage ;;
    *) usage ;;
  esac
done

[[ -n "${ANALYSIS}" && -n "${Q}" && -n "${S}" ]] || usage

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK=$(realpath -m "${WORK:-${ANALYSIS}/work}")
BLAST_DIR=$(realpath -m "${GENETRIBE_BLAST_DIR:-${WORK}/blast}")
DB_DIR="${BLAST_DIR}/db"
OUT="${BLAST_DIR}/${Q}_${S}.blast"
EVALUE="${GENETRIBE_EVALUE:-1e-5}"
THREADS="${GENETRIBE_THREADS:-${SLURM_CPUS_PER_TASK:-8}}"

long_q="${BLAST_DIR}/${Q}_long.fa"
db_s="${DB_DIR}/${S}"

if [[ "${GENETRIBE_SKIP_EXISTING:-1}" == "1" && -f "${OUT}" && -s "${OUT}" ]]; then
  # Treat empty as incomplete; non-empty skip
  echo "Skip existing ${Q} → ${S} (${OUT})"
  exit 0
fi

[[ -f "${long_q}" ]] || {
  echo "ERROR: missing ${long_q}; run: bash scripts/genetribe_prep_blast_dbs.sh ${ANALYSIS}" >&2
  exit 1
}
[[ -f "${db_s}.psq" || -f "${db_s}.pin" ]] || {
  echo "ERROR: missing BLAST DB ${db_s}; run: bash scripts/genetribe_prep_blast_dbs.sh ${ANALYSIS}" >&2
  exit 1
}
command -v blastp >/dev/null 2>&1 || {
  echo "ERROR: blastp not on PATH" >&2
  exit 1
}

mkdir -p "${BLAST_DIR}"
tmp="${OUT}.tmp.$$"
echo "=== blastp ${Q} → ${S} threads=${THREADS} evalue=${EVALUE} ==="
blastp \
  -query "${long_q}" \
  -db "${db_s}" \
  -evalue "${EVALUE}" \
  -num_threads "${THREADS}" \
  -outfmt 6 \
  -out "${tmp}"
mv "${tmp}" "${OUT}"
echo "Wrote ${OUT}"
