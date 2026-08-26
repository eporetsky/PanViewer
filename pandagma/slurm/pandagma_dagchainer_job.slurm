#!/bin/bash
# One standalone DAG batch job (not a job array). Submitted by submit_pandagma_dag_separate_jobs.sh
#
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION

set -euo pipefail

if [[ -z "${BATCH_MANIFEST:-}" || -z "${WORK:-}" ]]; then
  echo "ERROR: submit via scripts/submit_pandagma_dag_separate_jobs.sh (sets BATCH_MANIFEST, WORK)." >&2
  exit 1
fi

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then cd "${SLURM_SUBMIT_DIR}" || exit 1; fi

ENV_FILE="${WORK}/slurm/dag-job-env.sh"
[[ -f "${ENV_FILE}" ]] || { echo "ERROR: missing ${ENV_FILE}" >&2; exit 1; }
# shellcheck source=/dev/null
source "${ENV_FILE}"

export PATH="${REPO_ROOT}/pandagma-upstream/bin:${PATH}"
[[ -n "${PANDAGMA_CONDA_ENV:-}" ]] && export PATH="${PANDAGMA_CONDA_ENV}/bin:${PATH}"

echo "Job ${SLURM_JOB_ID} batch=${BATCH_ID:-?} manifest=${BATCH_MANIFEST}"
bash "${REPO_ROOT}/scripts/pandagma_dagchainer_batch_shard.sh" \
  -c "${CONF}" -d "${ANALYSIS}" -w "${WORK}" -m "${BATCH_MANIFEST}"
