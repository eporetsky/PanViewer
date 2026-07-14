#!/bin/bash
# One-shot ingest: pandagma pan -s ingest
# Usage: sbatch slurm/pandagma_ingest.slurm wheat

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --job-name=pdg-ingest
#SBATCH -N1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

set -euo pipefail

GENOME_ARG="${1:?usage: sbatch slurm/pandagma_ingest.slurm <genome e.g. wheat>}"
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then cd "${SLURM_SUBMIT_DIR}" || exit 1; fi
REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"

if [[ "${GENOME_ARG}" = /* ]]; then
  ANALYSIS=$(realpath "${GENOME_ARG}")
else
  ANALYSIS=$(realpath "${REPO_ROOT}/${GENOME_ARG}")
fi

_name=$(basename "${ANALYSIS}")
if [[ -n "${CONF:-}" ]]; then
  CONF=$(realpath "${CONF}")
elif [[ -f "${ANALYSIS}/config/pan.${_name}.conf" ]]; then
  CONF=$(realpath "${ANALYSIS}/config/pan.${_name}.conf")
elif [[ -f "${ANALYSIS}/config/${_name}_pan.conf" ]]; then
  CONF=$(realpath "${ANALYSIS}/config/${_name}_pan.conf")
else
  echo "ERROR: missing ${ANALYSIS}/config/pan.${_name}.conf or ${_name}_pan.conf (or set CONF=)." >&2
  exit 1
fi

WORK=$(realpath -m "${WORK:-${ANALYSIS}/work}")
mkdir -p "${WORK}" "${REPO_ROOT}/log"
export SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-16}"

export PATH="${REPO_ROOT}/pandagma-upstream/bin:${PATH}"
[[ -n "${PANDAGMA_CONDA_ENV:-}" ]] && export PATH="${PANDAGMA_CONDA_ENV}/bin:${PATH}"

pandagma pan -c "${CONF}" -d "${ANALYSIS}" -w "${WORK}" -n "${SLURM_CPUS_PER_TASK}" -s ingest
