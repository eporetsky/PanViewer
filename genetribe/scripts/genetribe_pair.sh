#!/usr/bin/env bash
# Run genetribe core for one accession pair.
#
# Usage:
#   genetribe_pair.sh -d <analysis_dir> -l <accA> -f <accB>
#
# Optional env:
#   WORK, GENETRIBE_THREADS, GENETRIBE_EVALUE, GENETRIBE_BSR,
#   GENETRIBE_SKIP_EXISTING, GENETRIBE_EXTRA, GENETRIBE_SKIP_ID_CHECK
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

WORK="${WORK:-${ANALYSIS}/work}"
mkdir -p "${WORK}"
WORK=$(cd "${WORK}" && pwd)
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
  [[ -f "${need}" ]] || {
    echo "ERROR: missing ${need}" >&2
    exit 1
  }
done

if [[ "${GENETRIBE_SKIP_EXISTING:-0}" == "1" && -f "${DONE}" ]]; then
  echo "Skip existing pair ${L} x ${F} (${DONE})"
  exit 0
fi

command -v genetribe >/dev/null 2>&1 || {
  if [[ -e "${REPO_ROOT}/genetribe-upstream/genetribe" ]]; then
    export PATH="${REPO_ROOT}/genetribe-upstream:${PATH}"
  else
    echo "ERROR: genetribe not on PATH. Run: make setup && conda activate genetribe" >&2
    exit 1
  fi
}

CHECK="${REPO_ROOT}/scripts/genetribe_check_ids.sh"
bash "${CHECK}" "${ACC_L}/${L}.fa" "${ACC_L}/${L}.bed"
bash "${CHECK}" "${ACC_F}/${F}.fa" "${ACC_F}/${F}.bed"

if [[ ! "${GENETRIBE_EXTRA:-}" =~ (^|[[:space:]])-s[[:space:]] ]]; then
  GENETRIBE_EXTRA="-s @ ${GENETRIBE_EXTRA:-}"
  GENETRIBE_EXTRA="${GENETRIBE_EXTRA%% }"
fi

mkdir -p "${PAIR_DIR}"
for bed in "${ACC_L}/${L}.bed" "${ACC_F}/${F}.bed"; do
  ncols=$(awk 'NF && $1 !~ /^#/ { print NF; exit }' "${bed}" || true)
  if [[ "${ncols}" != "6" ]]; then
    echo "ERROR: ${bed} has ${ncols:-0} columns; GeneTribe needs 6" >&2
    echo "  (chrom start end gene_id score strand). Re-run scripts/genetribe_prep.sh" >&2
    exit 1
  fi
done

ln -sfn "${ACC_L}/${L}.fa" "${PAIR_DIR}/${L}.fa"
ln -sfn "${ACC_L}/${L}.bed" "${PAIR_DIR}/${L}.bed"
ln -sfn "${ACC_L}/${L}.chrlist" "${PAIR_DIR}/${L}.chrlist"
ln -sfn "${ACC_F}/${F}.fa" "${PAIR_DIR}/${F}.fa"
ln -sfn "${ACC_F}/${F}.bed" "${PAIR_DIR}/${F}.bed"
ln -sfn "${ACC_F}/${F}.chrlist" "${PAIR_DIR}/${F}.chrlist"
[[ -f "${ACC_L}/${L}.confidence" ]] && ln -sfn "${ACC_L}/${L}.confidence" "${PAIR_DIR}/${L}.confidence"
[[ -f "${ACC_F}/${F}.confidence" ]] && ln -sfn "${ACC_F}/${F}.confidence" "${PAIR_DIR}/${F}.confidence"

cd "${PAIR_DIR}"
echo "=== genetribe core -l ${L} -f ${F} -n ${THREADS} (cwd=${PAIR_DIR}) ==="
# shellcheck disable=SC2086
genetribe core -l "${L}" -f "${F}" -n "${THREADS}" -e "${EVALUE}" -b "${BSR}" ${GENETRIBE_EXTRA:-}

[[ -f "${L}_${F}.RBH" || -f "${F}_${L}.RBH" ]] || {
  echo "ERROR: genetribe finished but no *.RBH under ${PAIR_DIR}" >&2
  if [[ -d "${PAIR_DIR}/genetribe_output" ]]; then
    echo "  genetribe_output present; last files:" >&2
    ls -lt "${PAIR_DIR}/genetribe_output" 2>/dev/null | head -20 >&2 || true
  fi
  exit 1
}

date -u +"done %Y-%m-%dT%H:%M:%SZ" >"${DONE}"
echo "Pair complete: ${L} x ${F}"
