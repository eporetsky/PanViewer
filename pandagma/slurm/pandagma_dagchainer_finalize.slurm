#!/bin/bash
# Build 05_filtered_pairs.tsv after all DAG shard jobs complete.
# Usage: sbatch slurm/pandagma_dagchainer_finalize.slurm barley

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --job-name=pdg-dag-fin
#SBATCH -N1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

set -euo pipefail

GENOME_ARG="${1:?usage: sbatch slurm/pandagma_dagchainer_finalize.slurm <genome>}"
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
  echo "ERROR: missing pan config." >&2
  exit 1
fi

WORK=$(realpath -m "${WORK:-${ANALYSIS}/work}")
mkdir -p "${REPO_ROOT}/log"

export PATH="${REPO_ROOT}/pandagma-upstream/bin:${PATH}"
[[ -n "${PANDAGMA_CONDA_ENV:-}" ]] && export PATH="${PANDAGMA_CONDA_ENV}/bin:${PATH}"

bash "${REPO_ROOT}/scripts/pandagma_dagchainer_finalize.sh" -c "${CONF}" -w "${WORK}"
bash "${REPO_ROOT}/scripts/audit-dagchainer.sh" "${GENOME_ARG}"
