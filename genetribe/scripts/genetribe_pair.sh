#!/usr/bin/env bash
# Run genetribe core for one accession pair.
#
# Usage:
#   genetribe_pair.sh -d <analysis_dir> -l <accA> -f <accB>
#
# Optional env:
#   WORK=              work root (default <analysis>/work)
#   GENETRIBE_BLAST_DIR=  precomputed BLAST dir for genetribe core -d
#                         (default <work>/blast). If all 4 files exist for the
#                         pair, GeneTribe skips BLAST entirely.
#   GENETRIBE_REQUIRE_BLAST=1  fail if the 4 directed BLASTs are missing
#   GENETRIBE_THREADS= blast/thread count (default SLURM_CPUS_PER_TASK or 8)
#   GENETRIBE_EVALUE=  BLAST e-value (default 1e-5)
#   GENETRIBE_BSR=     -b threshold (default 75)
#   GENETRIBE_SKIP_EXISTING=1  skip if done marker exists
#   GENETRIBE_EXTRA=   extra args for genetribe core (default includes -s @)
#   GENETRIBE_SKIP_ID_CHECK=1  skip BED↔FASTA gene-ID match check
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
  if [[ -x "${REPO_ROOT}/genetribe-upstream/genetribe" ]]; then
    export PATH="${REPO_ROOT}/genetribe-upstream:${PATH}"
  elif [[ -f "${REPO_ROOT}/genetribe-upstream/genetribe" ]]; then
    export PATH="${REPO_ROOT}/genetribe-upstream:${PATH}"
  else
    echo "ERROR: genetribe not on PATH. Run: make setup && conda activate genetribe" >&2
    exit 1
  fi
}

# Fail early if BED gene IDs and FASTA headers do not match (gene-level IDs required).
CHECK="${REPO_ROOT}/scripts/genetribe_check_id_separator.sh"
chmod +x "${CHECK}" 2>/dev/null || true
bash "${CHECK}" "${ACC_L}/${L}.fa" "${ACC_L}/${L}.bed"
bash "${CHECK}" "${ACC_F}/${F}.fa" "${ACC_F}/${F}.bed"

# Inputs are gene-level IDs: do not use GeneTribe default -s '.' (strips real dots in IDs).
if [[ ! "${GENETRIBE_EXTRA:-}" =~ (^|[[:space:]])-s[[:space:]] ]]; then
  GENETRIBE_EXTRA="-s @ ${GENETRIBE_EXTRA:-}"
  GENETRIBE_EXTRA="${GENETRIBE_EXTRA%% }"
fi

BLAST_DIR=$(realpath -m "${GENETRIBE_BLAST_DIR:-${WORK}/blast}")
# GeneTribe -d expects: A_B.blast B_A.blast A_A.blast B_B.blast
need_blasts=(
  "${BLAST_DIR}/${L}_${F}.blast"
  "${BLAST_DIR}/${F}_${L}.blast"
  "${BLAST_DIR}/${L}_${L}.blast"
  "${BLAST_DIR}/${F}_${F}.blast"
)
missing_blasts=0
for bf in "${need_blasts[@]}"; do
  if [[ ! -f "${bf}" || ! -s "${bf}" ]]; then
    missing_blasts=$((missing_blasts + 1))
  fi
done
if [[ "${GENETRIBE_REQUIRE_BLAST:-0}" == "1" && "${missing_blasts}" -gt 0 ]]; then
  echo "ERROR: GENETRIBE_REQUIRE_BLAST=1 but ${missing_blasts}/4 BLAST files missing under ${BLAST_DIR}" >&2
  echo "  Run: bash slurm/submit_genetribe_blasts.sh $(basename "${ANALYSIS}")" >&2
  echo "  Then: bash scripts/genetribe_check_blasts.sh $(basename "${ANALYSIS}")" >&2
  exit 1
fi
if [[ "${missing_blasts}" -eq 0 ]]; then
  echo "Using precomputed BLASTs from ${BLAST_DIR} (genetribe core -d)"
elif [[ -d "${BLAST_DIR}" ]]; then
  echo "WARN: ${missing_blasts}/4 BLASTs missing under ${BLAST_DIR}; GeneTribe will compute them" >&2
fi

mkdir -p "${PAIR_DIR}"
# GeneTribe expects prefix.fa / prefix.bed / prefix.chrlist in cwd
# Beds must be 6-col; refuse 7-col Pandagma layout early.
for bed in "${ACC_L}/${L}.bed" "${ACC_F}/${F}.bed"; do
  ncols=$(awk 'NF && $1 !~ /^#/ { print NF; exit }' "${bed}" || true)
  if [[ "${ncols}" != "6" ]]; then
    echo "ERROR: ${bed} has ${ncols:-0} columns; GeneTribe needs 6 (chrom start end gene_id score strand)." >&2
    echo "  From Pandagma 7-col beds run: bash scripts/genetribe_fix_beds.sh ${ANALYSIS}" >&2
    exit 1
  fi
done

ln -sfn "${ACC_L}/${L}.fa" "${PAIR_DIR}/${L}.fa"
ln -sfn "${ACC_L}/${L}.bed" "${PAIR_DIR}/${L}.bed"
ln -sfn "${ACC_L}/${L}.chrlist" "${PAIR_DIR}/${L}.chrlist"
ln -sfn "${ACC_F}/${F}.fa" "${PAIR_DIR}/${F}.fa"
ln -sfn "${ACC_F}/${F}.bed" "${PAIR_DIR}/${F}.bed"
ln -sfn "${ACC_F}/${F}.chrlist" "${PAIR_DIR}/${F}.chrlist"

# Optional confidence files
[[ -f "${ACC_L}/${L}.confidence" ]] && ln -sfn "${ACC_L}/${L}.confidence" "${PAIR_DIR}/${L}.confidence"
[[ -f "${ACC_F}/${F}.confidence" ]] && ln -sfn "${ACC_F}/${F}.confidence" "${PAIR_DIR}/${F}.confidence"

cd "${PAIR_DIR}"
echo "=== genetribe core -l ${L} -f ${F} -n ${THREADS} -d ${BLAST_DIR} (cwd=${PAIR_DIR}) ==="
# shellcheck disable=SC2086
genetribe core -l "${L}" -f "${F}" -n "${THREADS}" -e "${EVALUE}" -b "${BSR}" -d "${BLAST_DIR}" ${GENETRIBE_EXTRA:-}

[[ -f "${L}_${F}.RBH" || -f "${F}_${L}.RBH" ]] || {
  echo "ERROR: genetribe finished but no *.RBH under ${PAIR_DIR}" >&2
  echo "  (leftover genetribe_output usually means core died mid-run, often during jcvi ortholog)" >&2
  if [[ -d "${PAIR_DIR}/genetribe_output" ]]; then
    echo "  genetribe_output present; last files:" >&2
    ls -lt "${PAIR_DIR}/genetribe_output" 2>/dev/null | head -20 >&2 || true
  fi
  exit 1
}

date -u +"done %Y-%m-%dT%H:%M:%SZ" >"${DONE}"
echo "Pair complete: ${L} x ${F}"
