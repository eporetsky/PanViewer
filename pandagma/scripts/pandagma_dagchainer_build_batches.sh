#!/usr/bin/env bash
# Split all 04_dag/*_matches.tsv into N batch manifest files (~equal file count per batch).
#
# Usage:
#   pandagma_dagchainer_build_batches.sh -w WORK_DIR -n 20
#
# Writes: WORK_DIR/slurm/dag-batches/batch_000.txt (paths relative to WORK_DIR)
#
set -euo pipefail

WORK_DIR=""
NUM_BATCHES=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -w) WORK_DIR=$2; shift 2 ;;
    -n) NUM_BATCHES=$2; shift 2 ;;
    -h) head -12 "$0" | tail -11; exit 0 ;;
    *) echo "Unknown: $1" >&2; exit 1 ;;
  esac
done

[[ -n "${WORK_DIR}" && -n "${NUM_BATCHES}" ]] || {
  echo "ERROR: -w WORK_DIR -n NUM_BATCHES required" >&2
  exit 1
}
[[ "${NUM_BATCHES}" =~ ^[0-9]+$ ]] && [[ "${NUM_BATCHES}" -gt 0 ]] || {
  echo "ERROR: -n must be a positive integer" >&2
  exit 1
}

WORK_DIR=$(realpath "${WORK_DIR}")
DAG_DIR="${WORK_DIR}/04_dag"
OUT_DIR="${WORK_DIR}/slurm/dag-batches"

[[ -d "${DAG_DIR}" ]] || { echo "ERROR: missing ${DAG_DIR}" >&2; exit 1; }
mkdir -p "${OUT_DIR}"

mapfile -t all < <(
  for f in "${DAG_DIR}"/*_matches.tsv; do
    [[ -f "${f}" ]] || continue
    basename "${f}"
  done | sort
)
n=${#all[@]}
[[ "${n}" -gt 0 ]] || { echo "ERROR: no *_matches.tsv in ${DAG_DIR}" >&2; exit 1; }

rm -f "${OUT_DIR}"/batch_*.txt
i=0
for f in "${all[@]}"; do
  b=$((i % NUM_BATCHES))
  printf '04_dag/%s\n' "${f}" >> "${OUT_DIR}/batch_$(printf '%03d' "${b}").txt"
  i=$((i + 1))
done

echo "Split ${n} match files into ${NUM_BATCHES} batch manifest(s) under ${OUT_DIR}/"
for mf in "${OUT_DIR}"/batch_*.txt; do
  [[ -f "${mf}" ]] || continue
  printf '  %s: %s files\n' "$(basename "${mf}")" "$(wc -l <"${mf}" | tr -d ' ')"
done
