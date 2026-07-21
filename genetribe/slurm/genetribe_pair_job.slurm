#!/bin/bash
# One GeneTribe core job for a single accession pair.
# Submitted by slurm/submit_genetribe_pairs.sh — do not sbatch this directly
# without setting ANALYSIS, PAIR_L, PAIR_F, WORK, REPO_ROOT.
#
# Default resources (overridden by sbatch flags from the submitter):
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH -t 72:00:00
#SBATCH -N1
#SBATCH --ntasks=1

set -euo pipefail

if [[ -z "${ANALYSIS:-}" || -z "${PAIR_L:-}" || -z "${PAIR_F:-}" ]]; then
  echo "ERROR: submit via slurm/submit_genetribe_pairs.sh (sets ANALYSIS, PAIR_L, PAIR_F)." >&2
  exit 1
fi

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then cd "${SLURM_SUBMIT_DIR}" || exit 1; fi
REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"

# shellcheck source=/dev/null
source "${REPO_ROOT}/scripts/genetribe_env.sh"
_genetribe_bootstrap

export WORK="${WORK:-${ANALYSIS}/work}"
export GENETRIBE_BLAST_DIR="${GENETRIBE_BLAST_DIR:-${WORK}/blast}"
export SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-48}"
export GENETRIBE_THREADS="${GENETRIBE_THREADS:-${SLURM_CPUS_PER_TASK}}"
export GENETRIBE_SKIP_EXISTING="${GENETRIBE_SKIP_EXISTING:-1}"
export GENETRIBE_REQUIRE_BLAST="${GENETRIBE_REQUIRE_BLAST:-0}"

PAIR_DIR="${WORK}/pairs/${PAIR_L}_x_${PAIR_F}"
mkdir -p "${PAIR_DIR}" "${REPO_ROOT}/log"
PAIR_LOG="${PAIR_DIR}/pair.log"

{
  echo "=== GeneTribe pair job start $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
  echo "SLURM_JOB_ID=${SLURM_JOB_ID:-na} SLURM_JOB_NODELIST=${SLURM_JOB_NODELIST:-na}"
  echo "ANALYSIS=${ANALYSIS}"
  echo "PAIR=${PAIR_L} x ${PAIR_F}"
  echo "WORK=${WORK}"
  echo "BLAST_DIR=${GENETRIBE_BLAST_DIR} REQUIRE_BLAST=${GENETRIBE_REQUIRE_BLAST}"
  echo "PAIR_DIR=${PAIR_DIR}"
  echo "CPUS=${SLURM_CPUS_PER_TASK} GENETRIBE_THREADS=${GENETRIBE_THREADS}"
  echo "MEM_HINT=${SLURM_MEM_PER_NODE:-${SLURM_MEM_PER_CPU:-na}} SKIP_EXISTING=${GENETRIBE_SKIP_EXISTING}"
  echo "which genetribe: $(command -v genetribe || echo MISSING)"
} | tee -a "${PAIR_LOG}"

bash "${REPO_ROOT}/scripts/genetribe_pair.sh" -d "${ANALYSIS}" -l "${PAIR_L}" -f "${PAIR_F}" 2>&1 | tee -a "${PAIR_LOG}"
rc=${PIPESTATUS[0]}

{
  echo "=== GeneTribe pair job end $(date -u +%Y-%m-%dT%H:%M:%SZ) exit=${rc} ==="
} | tee -a "${PAIR_LOG}"

exit "${rc}"
