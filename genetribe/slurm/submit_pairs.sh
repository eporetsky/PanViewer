#!/usr/bin/env bash
# Submit one Slurm job per unordered accession pair (GeneTribe core).
#
# Usage (from genetribe/):
#   bash slurm/submit_pairs.sh wheat
#
# Override with env: GT_CPUS, GT_MEM, GT_TIME, SLURM_ACCOUNT, SLURM_PARTITION, …
set -euo pipefail

GENOME_ARG="${1:?usage: $0 <genome e.g. wheat>}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ "${GENOME_ARG}" = /* ]]; then
  ANALYSIS=$(realpath "${GENOME_ARG}")
else
  ANALYSIS=$(realpath "${REPO_ROOT}/${GENOME_ARG}")
fi

WORK="${WORK:-${ANALYSIS}/work}"
mkdir -p "${REPO_ROOT}/log" "${WORK}/slurm" "${WORK}/pairs"
WORK=$(cd "${WORK}" && pwd)

SUBMIT_LOG="${REPO_ROOT}/log/submit_pairs.$(date -u +%Y%m%dT%H%M%SZ).log"
exec > >(tee -a "${SUBMIT_LOG}") 2>&1
echo "=== submit_pairs $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
echo "REPO_ROOT=${REPO_ROOT}"
echo "ANALYSIS=${ANALYSIS}"
echo "WORK=${WORK}"
echo "submit log: ${SUBMIT_LOG}"

bash "${REPO_ROOT}/scripts/genetribe_check_ids.sh" "${ANALYSIS}"

if [[ ! "${GENETRIBE_EXTRA:-}" =~ (^|[[:space:]])-s[[:space:]] ]]; then
  GENETRIBE_EXTRA="-s @ ${GENETRIBE_EXTRA:-}"
  GENETRIBE_EXTRA="${GENETRIBE_EXTRA%% }"
fi

PAIRS_TSV="${WORK}/slurm/pairs.tsv"
bash "${REPO_ROOT}/scripts/genetribe_build_pairs.sh" "${ANALYSIS}" >"${PAIRS_TSV}"
n_pairs=$(wc -l <"${PAIRS_TSV}" | tr -d ' ')
[[ "${n_pairs}" -gt 0 ]] || { echo "ERROR: no pairs generated" >&2; exit 1; }

GT_CPUS="${GT_CPUS:-48}"
GT_MEM="${GT_MEM:-128G}"
GT_TIME="${GT_TIME:-72:00:00}"
GENETRIBE_SKIP_EXISTING="${GENETRIBE_SKIP_EXISTING:-1}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-small_grains}"
SLURM_PARTITION="${SLURM_PARTITION:-atlas}"
: "${GENETRIBE_HOME:=${REPO_ROOT}/genetribe-upstream}"

if [[ ! -e "${GENETRIBE_HOME}/genetribe" ]]; then
  echo "ERROR: ${GENETRIBE_HOME}/genetribe missing. From genetribe/ run: make setup" >&2
  exit 1
fi

JOB_SLURM="${REPO_ROOT}/slurm/pair_job.slurm"
echo "=== ${n_pairs} GeneTribe pairs (SKIP_EXISTING=${GENETRIBE_SKIP_EXISTING}) ==="
echo "  account=${SLURM_ACCOUNT} partition=${SLURM_PARTITION}"
echo "  per job: CPUS=${GT_CPUS} MEM=${GT_MEM} TIME=${GT_TIME}"

submitted=0
skipped=0
job_ids=()
while IFS=$'\t' read -r a b; do
  [[ -z "${a:-}" || -z "${b:-}" ]] && continue
  safe="${a}_x_${b}"
  safe="${safe//[^A-Za-z0-9._-]/_}"
  pair_dir="${WORK}/pairs/${a}_x_${b}"

  if [[ "${GENETRIBE_SKIP_EXISTING}" == "1" ]]; then
    if [[ -f "${pair_dir}/.done" ]] \
      || [[ -f "${pair_dir}/${a}_${b}.RBH" ]] \
      || [[ -f "${pair_dir}/${b}_${a}.RBH" ]]; then
      echo "  skip ${a} x ${b} (already complete)"
      skipped=$((skipped + 1))
      continue
    fi
  fi

  jid=$(sbatch --parsable \
    --account="${SLURM_ACCOUNT}" \
    --partition="${SLURM_PARTITION}" \
    --job-name="gt-${safe}" \
    --cpus-per-task="${GT_CPUS}" \
    --mem="${GT_MEM}" \
    -t "${GT_TIME}" \
    -N1 \
    --export=ALL,ANALYSIS="${ANALYSIS}",WORK="${WORK}",REPO_ROOT="${REPO_ROOT}",PAIR_L="${a}",PAIR_F="${b}",GENETRIBE_SKIP_EXISTING="${GENETRIBE_SKIP_EXISTING}",GENETRIBE_THREADS="${GT_CPUS}",GENETRIBE_CONDA_ENV="${GENETRIBE_CONDA_ENV:-}",GENETRIBE_HOME="${GENETRIBE_HOME}",GENETRIBE_EXTRA="${GENETRIBE_EXTRA:-}" \
    -o "${REPO_ROOT}/log/stdout.${safe}.%j.%N" \
    -e "${REPO_ROOT}/log/stderr.${safe}.%j.%N" \
    "${JOB_SLURM}")
  echo "  ${a} x ${b}: job ${jid}"
  job_ids+=("${jid}")
  submitted=$((submitted + 1))
done <"${PAIRS_TSV}"

echo "${job_ids[*]}" >"${WORK}/slurm/submitted_job_ids.txt"
echo ""
echo "Submitted ${submitted} job(s); skipped ${skipped} already-complete pair(s)."
echo "When all finish:"
echo "  sbatch slurm/cluster_pans.slurm ${GENOME_ARG}"
echo "  sbatch slurm/build_db.slurm ${GENOME_ARG}"
