#!/bin/bash
# After ingest, sharded mmseqs+filter, and sharded DAG (see pandagma_dagchainer_array.slurm):
# run pandagma pan from mcl through tabularize.
#
# Full barley DAG + pan workflow:
#   bash slurm/submit_pandagma_dag_separate_jobs.sh barley
#   # wait for all pdg-dag jobs
#   sbatch slurm/pandagma_dagchainer_finalize.slurm barley
#   sbatch slurm/pandagma_pan_resume.slurm barley
#
# Usage: sbatch slurm/pandagma_pan_resume.slurm barley
# Do NOT run pandagma pan -s dagchainer if you used the DAG array (it redoes all pairs on one node).

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --job-name=pdg-pan-resume
#SBATCH -N1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=256GB
#SBATCH -t 5-00:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

set -euo pipefail

GENOME_ARG="${1:?usage: sbatch slurm/pandagma_pan_resume.slurm <genome e.g. barley>}"
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
export SLURM_CPUS_PER_TASK="${SLURM_CPUS_PER_TASK:-48}"

export PATH="${REPO_ROOT}/pandagma-upstream/bin:${PATH}"
[[ -n "${PANDAGMA_CONDA_ENV:-}" ]] && export PATH="${PANDAGMA_CONDA_ENV}/bin:${PATH}"

if [[ ! -f "${WORK}/05_filtered_pairs.tsv" ]]; then
  echo "========== build 05_filtered_pairs.tsv (DAG finalize) =========="
  bash "${REPO_ROOT}/scripts/pandagma_dagchainer_finalize.sh" -c "${CONF}" -w "${WORK}"
fi

steps=(
  mcl consense cluster_rest add_extra
  check_leftovers tabularize
)

for s in "${steps[@]}"; do
  echo "========== pandagma pan -s ${s} =========="
  pandagma pan -c "${CONF}" -d "${ANALYSIS}" -w "${WORK}" -n "${SLURM_CPUS_PER_TASK}" -s "${s}"
done

echo "Resume done through tabularize. conf=${CONF} work=${WORK}"
