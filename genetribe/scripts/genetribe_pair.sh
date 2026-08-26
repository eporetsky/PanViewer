#!/usr/bin/env bash
# Run genetribe core for one accession pair.
#
# Usage:
#   genetribe_pair.sh -d <analysis_dir> -l <accA> -f <accB>
#
# Optional env:
#   GENETRIBE_BLAST_DIR=  precomputed BLAST directory (default <work>/blast)
#   GENETRIBE_THREADS=    default SLURM_CPUS_PER_TASK or 8
#   GENETRIBE_SKIP_EXISTING=1
set -euo pipefail

ANALYSIS=""
L=""
F=""

usage() {
  echo "usage: $0 -d <analysis_dir> -l <accA> -f <accB>" >&2
  exit 1
}

while getopts ":d:l:f:h" opt; do
  case "${opt}" in
    d) ANALYSIS=$(realpath "${OPTARG}") ;;
    l) L="${OPTARG}" ;;
    f) F="${OPTARG}" ;;
    h) usage ;;
    *) usage ;;
  esac
done

[[ -n "${ANALYSIS}" && -n "${L}" && -n "${F}" ]] || usage

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK=$(realpath -m "${WORK:-${ANALYSIS}/work}")
PAIR_DIR="${WORK}/pairs/${L}_x_${F}"
DONE="${PAIR_DIR}/.done"
THREADS="${GENETRIBE_THREADS:-${SLURM_CPUS_PER_TASK:-8}}"
EVALUE="${GENETRIBE_EVALUE:-1e-5}"
BSR="${GENETRIBE_BSR:-75}"

ACC_L="${ANALYSIS}/accessions/${L}"
ACC_F="${ANALYSIS}/accessions/${F}"
for need in \
  "${ACC_L}/${L}.fa" "${ACC_L}/${L}.bed" "${ACC_L}/${L}.chrlist" \
  "${ACC_F}/${F}.fa" "${ACC_F}/${F}.bed" "${ACC_F}/${F}.chrlist"
do
  [[ -f "${need}" ]] || { echo "ERROR: missing ${need}" >&2; exit 1; }
done

if [[ "${GENETRIBE_SKIP_EXISTING:-0}" == "1" && -f "${DONE}" ]]; then
  echo "Skip existing pair ${L} x ${F}"
  exit 0
fi

# shellcheck source=/dev/null
source "${REPO_ROOT}/scripts/genetribe_env.sh"
_genetribe_bootstrap

bash "${REPO_ROOT}/scripts/genetribe_check_ids.sh" "${ACC_L}/${L}.fa" "${ACC_L}/${L}.bed"
bash "${REPO_ROOT}/scripts/genetribe_check_ids.sh" "${ACC_F}/${F}.fa" "${ACC_F}/${F}.bed"

if [[ ! "${GENETRIBE_EXTRA:-}" =~ (^|[[:space:]])-s[[:space:]] ]]; then
  GENETRIBE_EXTRA="-s @ ${GENETRIBE_EXTRA:-}"
  GENETRIBE_EXTRA="${GENETRIBE_EXTRA%% }"
fi

BLAST_DIR=$(realpath -m "${GENETRIBE_BLAST_DIR:-${WORK}/blast}")
mkdir -p "${PAIR_DIR}"

for bed in "${ACC_L}/${L}.bed" "${ACC_F}/${F}.bed"; do
  ncols=$(awk 'NF && $1 !~ /^#/ { print NF; exit }' "${bed}" || true)
  if [[ "${ncols}" != "6" ]]; then
    echo "ERROR: ${bed} has ${ncols:-0} columns; GeneTribe needs 6" >&2
    exit 1
  fi
done

ln -sfn "${ACC_L}/${L}.fa" "${PAIR_DIR}/${L}.fa"
ln -sfn "${ACC_L}/${L}.bed" "${PAIR_DIR}/${L}.bed"
ln -sfn "${ACC_L}/${L}.chrlist" "${PAIR_DIR}/${L}.chrlist"
ln -sfn "${ACC_F}/${F}.fa" "${PAIR_DIR}/${F}.fa"
ln -sfn "${ACC_F}/${F}.bed" "${PAIR_DIR}/${F}.bed"
ln -sfn "${ACC_F}/${F}.chrlist" "${PAIR_DIR}/${F}.chrlist"

cd "${PAIR_DIR}"
echo "=== genetribe core -l ${L} -f ${F} -n ${THREADS} -d ${BLAST_DIR} ==="
# shellcheck disable=SC2086
genetribe core -l "${L}" -f "${F}" -n "${THREADS}" -e "${EVALUE}" -b "${BSR}" -d "${BLAST_DIR}" ${GENETRIBE_EXTRA:-}

[[ -f "${L}_${F}.RBH" || -f "${F}_${L}.RBH" ]] || {
  echo "ERROR: no *.RBH under ${PAIR_DIR}" >&2
  exit 1
}

date -u +"done %Y-%m-%dT%H:%M:%SZ" >"${DONE}"
echo "Pair complete: ${L} x ${F}"
