#!/bin/bash
# One directed BLASTP job (QUERY → SUBJECT) for GeneTribe -d reuse.
# Submitted by slurm/submit_genetribe_blasts.sh.
#
# Default resources (overridden by sbatch flags from the submitter):
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -N1
#SBATCH --ntasks=1

set -euo pipefail

if [[ -z "${ANALYSIS:-}" || -z "${BLAST_Q:-}" || -z "${BLAST_S:-}" ]]; then
  echo "ERROR: submit via slurm/submit_genetribe_blasts.sh (sets ANALYSIS, BLAST_Q, BLAST_S)." >&2
  exit 1
fi

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then cd "${SLURM_SUBMIT_DIR}" || exit 1; fi
REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"

# shellcheck source=/dev/null
source "${REPO_ROOT}/scripts/genetribe_env.sh"
_genetribe_bootstrap

export WORK="${WORK:-${ANALYSIS}/work}"
export GENETRIBE_BLAST_DIR="${GENETRIBE_BLAST_DIR:-${WORK}/blast}"
export SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-24}"
export GENETRIBE_THREADS="${GENETRIBE_THREADS:-${SLURM_CPUS_PER_TASK}}"
export GENETRIBE_SKIP_EXISTING="${GENETRIBE_SKIP_EXISTING:-1}"

echo "=== BLAST ${BLAST_Q} → ${BLAST_S} start $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
echo "ANALYSIS=${ANALYSIS} BLAST_DIR=${GENETRIBE_BLAST_DIR} THREADS=${GENETRIBE_THREADS}"

bash "${REPO_ROOT}/scripts/genetribe_blast_one.sh" \
  -d "${ANALYSIS}" -q "${BLAST_Q}" -s "${BLAST_S}"
rc=$?

echo "=== BLAST ${BLAST_Q} → ${BLAST_S} end $(date -u +%Y-%m-%dT%H:%M:%SZ) exit=${rc} ==="
exit "${rc}"
