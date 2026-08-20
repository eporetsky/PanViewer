#!/bin/bash
# Cluster pan-genes: subgenome-aware + collinear RBH.
#
# Usage:
#   sbatch slurm/cluster_pans.slurm wheat
#
# Writes:
#   <species>/work/genetribe_subgenome_pans.{hsh,clust}.tsv
#
# Then:
#   sbatch slurm/build_db.slurm <species>
#
#SBATCH --account=small_grains
#SBATCH --partition=atlas
#SBATCH --job-name=gt-cluster
#SBATCH -N1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH -t 08:00:00
#SBATCH -o "./log/stdout.cluster_pans.%j.%N"
#SBATCH -e "./log/stderr.cluster_pans.%j.%N"

set -euo pipefail

GENOME_ARG="${1:?usage: sbatch slurm/cluster_pans.slurm <species>}"
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then cd "${SLURM_SUBMIT_DIR}" || exit 1; fi
REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
mkdir -p "${REPO_ROOT}/log"


if [[ "${GENOME_ARG}" = /* ]]; then
  ANALYSIS=$(realpath "${GENOME_ARG}")
else
  ANALYSIS=$(realpath "${REPO_ROOT}/${GENOME_ARG}")
fi

# Optional account/partition from env are set in #SBATCH headers above;
# override at submit time: sbatch --account=… --partition=… slurm/cluster_pans.slurm …

# shellcheck source=/dev/null
source "${REPO_ROOT}/scripts/genetribe_env.sh"
_genetribe_bootstrap || true

LOG="${ANALYSIS}/work/cluster_pans.log"
{
  echo "=== cluster pans $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
  echo "ANALYSIS=${ANALYSIS}"
  echo "SLURM_JOB_ID=${SLURM_JOB_ID:-na} node=${SLURM_NODELIST:-na}"
  python3 "${REPO_ROOT}/scripts/genetribe_cluster_pans.py" "${ANALYSIS}"
  echo "Done $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} 2>&1 | tee "${LOG}"
