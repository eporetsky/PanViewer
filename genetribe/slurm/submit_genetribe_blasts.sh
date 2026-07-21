#!/usr/bin/env bash
# Precompute all directed BLASTPs (including self) for GeneTribe core -d reuse.
#
# Usage:
#   bash slurm/submit_genetribe_blasts.sh <analysis_dir>
#
# After BLAST jobs finish:
#   bash scripts/genetribe_check_blasts.sh <analysis_dir>
#   bash slurm/submit_genetribe_pairs.sh <analysis_dir>
#
# Resources / site (override as needed):
#   GT_CPUS=24 GT_MEM=64G GT_TIME=12:00:00
#   GT_ACCOUNT=myaccount GT_PARTITION=mypartition
#   GENETRIBE_CONDA_ENV=/path/to/conda/envs/genetribe
#   GENETRIBE_SKIP_EXISTING=1   # default: skip non-empty *.blast
#
set -euo pipefail

GENOME_ARG="${1:?usage: $0 <analysis_dir>}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ "${GENOME_ARG}" = /* ]]; then
  ANALYSIS=$(realpath "${GENOME_ARG}")
else
  ANALYSIS=$(realpath "${REPO_ROOT}/${GENOME_ARG}")
fi

WORK=$(realpath -m "${WORK:-${ANALYSIS}/work}")
BLAST_DIR=$(realpath -m "${GENETRIBE_BLAST_DIR:-${WORK}/blast}")
mkdir -p "${REPO_ROOT}/log" "${WORK}/slurm" "${BLAST_DIR}"

SUBMIT_LOG="${REPO_ROOT}/log/submit_blasts.$(date -u +%Y%m%dT%H%M%SZ).log"
exec > >(tee -a "${SUBMIT_LOG}") 2>&1
echo "=== submit_genetribe_blasts $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
echo "REPO_ROOT=${REPO_ROOT}"
echo "ANALYSIS=${ANALYSIS}"
echo "WORK=${WORK}"
echo "BLAST_DIR=${BLAST_DIR}"
echo "submit log: ${SUBMIT_LOG}"

: "${GENETRIBE_HOME:=${REPO_ROOT}/genetribe-upstream}"
export REPO_ROOT GENETRIBE_HOME
export GENETRIBE_CONDA_ENV="${GENETRIBE_CONDA_ENV:-}"
export GENETRIBE_EXTRA="${GENETRIBE_EXTRA:-}"

# Ensure -s @ unless user set -s (must match pair jobs / longestfasta)
if [[ ! "${GENETRIBE_EXTRA:-}" =~ (^|[[:space:]])-s[[:space:]] ]]; then
  GENETRIBE_EXTRA="-s @ ${GENETRIBE_EXTRA:-}"
  GENETRIBE_EXTRA="${GENETRIBE_EXTRA%% }"
fi
export GENETRIBE_EXTRA

# shellcheck source=/dev/null
source "${REPO_ROOT}/scripts/genetribe_env.sh"
_genetribe_bootstrap
_genetribe_sbatch_site_args

BUILD="${REPO_ROOT}/scripts/genetribe_build_blast_tasks.sh"
PREP_DB="${REPO_ROOT}/scripts/genetribe_prep_blast_dbs.sh"
ID_CHECK="${REPO_ROOT}/scripts/genetribe_check_id_separator.sh"
chmod +x \
  "${BUILD}" "${PREP_DB}" "${ID_CHECK}" \
  "${REPO_ROOT}/scripts/genetribe_blast_one.sh" \
  "${REPO_ROOT}/scripts/genetribe_blast_sep.sh" \
  "${REPO_ROOT}/scripts/genetribe_check_blasts.sh" \
  "${REPO_ROOT}/scripts/genetribe_list_accessions.sh" \
  2>/dev/null || true

bash "${ID_CHECK}" "${ANALYSIS}"

echo "=== preparing longest FASTA + BLAST DBs ==="
bash "${PREP_DB}" "${ANALYSIS}"

TASKS_TSV="${WORK}/slurm/blast_tasks.tsv"
bash "${BUILD}" "${ANALYSIS}" >"${TASKS_TSV}"
n_tasks=$(wc -l <"${TASKS_TSV}" | tr -d ' ')
[[ "${n_tasks}" -gt 0 ]] || {
  echo "ERROR: no BLAST tasks generated" >&2
  exit 1
}

GT_CPUS="${GT_CPUS:-24}"
GT_MEM="${GT_MEM:-64G}"
GT_TIME="${GT_TIME:-12:00:00}"
GENETRIBE_SKIP_EXISTING="${GENETRIBE_SKIP_EXISTING:-1}"

JOB_SLURM="${REPO_ROOT}/slurm/genetribe_blast_job.slurm"
echo "=== ${n_tasks} directed BLAST tasks (SKIP_EXISTING=${GENETRIBE_SKIP_EXISTING}) ==="
echo "  account=${GT_ACCOUNT:-<default>} partition=${GT_PARTITION:-<default>}"
echo "  per job: CPUS=${GT_CPUS} MEM=${GT_MEM} TIME=${GT_TIME}"

submitted=0
skipped=0
job_ids=()
while IFS=$'\t' read -r q s; do
  [[ -z "${q:-}" || -z "${s:-}" ]] && continue
  out="${BLAST_DIR}/${q}_${s}.blast"
  if [[ "${GENETRIBE_SKIP_EXISTING}" == "1" && -f "${out}" && -s "${out}" ]]; then
    skipped=$((skipped + 1))
    continue
  fi

  safe="${q}_v_${s}"
  safe="${safe//[^A-Za-z0-9._-]/_}"
  jid=$(sbatch --parsable \
    "${SBATCH_SITE_ARGS[@]}" \
    --job-name="gtb-${safe}" \
    --cpus-per-task="${GT_CPUS}" \
    --mem="${GT_MEM}" \
    -t "${GT_TIME}" \
    -N1 \
    --export=ALL,ANALYSIS="${ANALYSIS}",WORK="${WORK}",REPO_ROOT="${REPO_ROOT}",BLAST_Q="${q}",BLAST_S="${s}",GENETRIBE_BLAST_DIR="${BLAST_DIR}",GENETRIBE_SKIP_EXISTING="${GENETRIBE_SKIP_EXISTING}",GENETRIBE_THREADS="${GT_CPUS}",GENETRIBE_CONDA_ENV="${GENETRIBE_CONDA_ENV:-}",GENETRIBE_HOME="${GENETRIBE_HOME}",GENETRIBE_EXTRA="${GENETRIBE_EXTRA:-}" \
    -o "${REPO_ROOT}/log/stdout.blast.${safe}.%j.%N" \
    -e "${REPO_ROOT}/log/stderr.blast.${safe}.%j.%N" \
    "${JOB_SLURM}")
  echo "  ${q} → ${s}: job ${jid}"
  job_ids+=("${jid}")
  submitted=$((submitted + 1))
done <"${TASKS_TSV}"

echo "${job_ids[*]}" >"${WORK}/slurm/submitted_blast_job_ids.txt"
echo ""
echo "Submitted ${submitted} BLAST job(s); skipped ${skipped} existing."
echo "  IDs: ${WORK}/slurm/submitted_blast_job_ids.txt"
echo "Queue: squeue -u \$USER"
echo "When all finish:"
echo "  bash scripts/genetribe_check_blasts.sh ${GENOME_ARG}"
echo "  bash slurm/submit_genetribe_pairs.sh ${GENOME_ARG}"
