#!/usr/bin/env bash
# Run all unordered accession pairs locally (no Slurm).
# For clusters, prefer: bash slurm/submit_pairs.sh <species>
#
# Usage:
#   bash scripts/genetribe_run_pairs.sh wheat
#   GENETRIBE_THREADS=16 bash scripts/genetribe_run_pairs.sh wheat
#
# Optional: MAX_PARALLEL=4  (default 1 = serial)
set -euo pipefail

GENOME_ARG="${1:?usage: $0 <genome e.g. wheat>}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ "${GENOME_ARG}" = /* ]]; then
  ANALYSIS=$(realpath "${GENOME_ARG}")
else
  ANALYSIS=$(realpath "${REPO_ROOT}/${GENOME_ARG}")
fi

WORK="${WORK:-${ANALYSIS}/work}"
mkdir -p "${WORK}/slurm" "${WORK}/pairs" "${REPO_ROOT}/log"
WORK=$(cd "${WORK}" && pwd)

export REPO_ROOT
# shellcheck source=/dev/null
source "${REPO_ROOT}/scripts/genetribe_env.sh"
_genetribe_bootstrap

bash "${REPO_ROOT}/scripts/genetribe_check_ids.sh" "${ANALYSIS}"

if [[ ! "${GENETRIBE_EXTRA:-}" =~ (^|[[:space:]])-s[[:space:]] ]]; then
  GENETRIBE_EXTRA="-s @ ${GENETRIBE_EXTRA:-}"
  GENETRIBE_EXTRA="${GENETRIBE_EXTRA%% }"
fi
export GENETRIBE_EXTRA
export GENETRIBE_SKIP_EXISTING="${GENETRIBE_SKIP_EXISTING:-1}"
export GENETRIBE_THREADS="${GENETRIBE_THREADS:-8}"

PAIRS_TSV="${WORK}/slurm/pairs.tsv"
bash "${REPO_ROOT}/scripts/genetribe_build_pairs.sh" "${ANALYSIS}" >"${PAIRS_TSV}"
n_pairs=$(wc -l <"${PAIRS_TSV}" | tr -d ' ')
[[ "${n_pairs}" -gt 0 ]] || { echo "ERROR: no pairs" >&2; exit 1; }

MAX_PARALLEL="${MAX_PARALLEL:-1}"
echo "=== ${n_pairs} pairs (local, MAX_PARALLEL=${MAX_PARALLEL}, threads=${GENETRIBE_THREADS}) ==="

run_one() {
  local a="$1" b="$2"
  local pair_dir="${WORK}/pairs/${a}_x_${b}"
  if [[ "${GENETRIBE_SKIP_EXISTING}" == "1" ]]; then
    if [[ -f "${pair_dir}/.done" ]] \
      || [[ -f "${pair_dir}/${a}_${b}.RBH" ]] \
      || [[ -f "${pair_dir}/${b}_${a}.RBH" ]]; then
      echo "  skip ${a} x ${b}"
      return 0
    fi
  fi
  echo "  run ${a} x ${b}"
  bash "${REPO_ROOT}/scripts/genetribe_pair.sh" -d "${ANALYSIS}" -l "${a}" -f "${b}"
}

if [[ "${MAX_PARALLEL}" -le 1 ]]; then
  while IFS=$'\t' read -r a b; do
    [[ -z "${a:-}" || -z "${b:-}" ]] && continue
    run_one "${a}" "${b}"
  done <"${PAIRS_TSV}"
else
  # Simple background pool
  running=0
  while IFS=$'\t' read -r a b; do
    [[ -z "${a:-}" || -z "${b:-}" ]] && continue
    while [[ "${running}" -ge "${MAX_PARALLEL}" ]]; do
      wait -n 2>/dev/null || wait
      running=$((running - 1))
    done
    run_one "${a}" "${b}" &
    running=$((running + 1))
  done <"${PAIRS_TSV}"
  wait
fi

echo "All pairs finished. Next:"
echo "  python3 scripts/genetribe_cluster_pans.py ${ANALYSIS}"
echo "  bash scripts/genetribe_build_db.sh ${GENOME_ARG}"
