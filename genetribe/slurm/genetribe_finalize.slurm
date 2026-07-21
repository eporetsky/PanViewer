#!/bin/bash
# After all GeneTribe pair jobs finish: RBH connected components → PanViewer TSV.
# Usage:
#   sbatch [--account=...] [--partition=...] slurm/genetribe_finalize.slurm <analysis_dir>
#
# Resources: light CPU; raise mem if finalize OOMs on large graphs.
#
#SBATCH --job-name=gt-finalize
#SBATCH -N1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o "./log/stdout.finalize.%j.%N"
#SBATCH -e "./log/stderr.finalize.%j.%N"

set -euo pipefail

GENOME_ARG="${1:?usage: sbatch slurm/genetribe_finalize.slurm <analysis_dir>}"
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then cd "${SLURM_SUBMIT_DIR}" || exit 1; fi
REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
mkdir -p "${REPO_ROOT}/log"

if [[ "${GENOME_ARG}" = /* ]]; then
  ANALYSIS=$(realpath "${GENOME_ARG}")
else
  ANALYSIS=$(realpath "${REPO_ROOT}/${GENOME_ARG}")
fi

WORK=$(realpath -m "${WORK:-${ANALYSIS}/work}")
mkdir -p "${WORK}"
FIN_LOG="${WORK}/finalize.log"

# shellcheck source=/dev/null
source "${REPO_ROOT}/scripts/genetribe_env.sh"
_genetribe_bootstrap || true   # finalize only needs python3; genetribe optional

{
  echo "=== GeneTribe finalize start $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
  echo "SLURM_JOB_ID=${SLURM_JOB_ID:-na} node=${SLURM_JOB_NODELIST:-na}"
  echo "ANALYSIS=${ANALYSIS} WORK=${WORK}"
  echo "CPUS=${SLURM_CPUS_PER_TASK:-4}"
} | tee -a "${FIN_LOG}"

PREFIX="${OGI_PREFIX:-OGI}"
python3 "${REPO_ROOT}/scripts/genetribe_rbh_to_pans.py" "${ANALYSIS}" --work "${WORK}" --prefix "${PREFIX}" 2>&1 | tee -a "${FIN_LOG}"
rc=${PIPESTATUS[0]}

{
  echo "Finalize exit=${rc} $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "PanViewer inputs:"
  echo "  ${WORK}/genetribe_pans.hsh.tsv"
  echo "  ${WORK}/genetribe_pans.clust.tsv"
  echo "Copy the hsh (or clust) TSV plus bed/ into PanViewer input/<species>/ then:"
  echo "  python build_index.py --force"
} | tee -a "${FIN_LOG}"

exit "${rc}"
