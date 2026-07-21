#!/usr/bin/env bash
# Submit one Slurm job per unordered accession pair (GeneTribe core).
#
# Usage:
#   bash slurm/submit_genetribe_pairs.sh <analysis_dir>
#
# Prefer precomputed BLASTs first:
#   bash slurm/submit_genetribe_blasts.sh <analysis_dir>
#   bash scripts/genetribe_check_blasts.sh <analysis_dir>
#   bash slurm/submit_genetribe_pairs.sh <analysis_dir>
#
# Resources / site (override as needed):
#   GT_CPUS GT_MEM GT_TIME
#   GT_ACCOUNT GT_PARTITION
#   GENETRIBE_CONDA_ENV=/path/to/conda/envs/genetribe
#   GENETRIBE_SKIP_EXISTING=1   # skip pairs with .done or *.RBH
#   GENETRIBE_EXTRA='-s @'      # applied by default if -s unset
#
# After all pairs finish:
#   sbatch [--account=...] [--partition=...] slurm/genetribe_finalize.slurm <analysis_dir>
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
mkdir -p "${REPO_ROOT}/log" "${WORK}/slurm" "${WORK}/pairs"

SUBMIT_LOG="${REPO_ROOT}/log/submit_pairs.$(date -u +%Y%m%dT%H%M%SZ).log"
exec > >(tee -a "${SUBMIT_LOG}") 2>&1
echo "=== submit_genetribe_pairs $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
echo "REPO_ROOT=${REPO_ROOT}"
echo "ANALYSIS=${ANALYSIS}"
echo "WORK=${WORK}"
echo "submit log: ${SUBMIT_LOG}"

BUILD="${REPO_ROOT}/scripts/genetribe_build_pairs.sh"
PAIR_RUNNER="${REPO_ROOT}/scripts/genetribe_pair.sh"
ID_CHECK="${REPO_ROOT}/scripts/genetribe_check_id_separator.sh"
chmod +x "${BUILD}" "${PAIR_RUNNER}" "${ID_CHECK}" "${REPO_ROOT}/scripts/genetribe_list_accessions.sh" 2>/dev/null || true

# shellcheck source=/dev/null
source "${REPO_ROOT}/scripts/genetribe_env.sh"
_genetribe_sbatch_site_args

# Preflight: gene-level BED↔FASTA IDs (no GeneTribe transcript stripping).
bash "${ID_CHECK}" "${ANALYSIS}"

# Default -s @ for all submit jobs unless user already set -s in GENETRIBE_EXTRA
if [[ ! "${GENETRIBE_EXTRA:-}" =~ (^|[[:space:]])-s[[:space:]] ]]; then
  GENETRIBE_EXTRA="-s @ ${GENETRIBE_EXTRA:-}"
  GENETRIBE_EXTRA="${GENETRIBE_EXTRA%% }"
fi

PAIRS_TSV="${WORK}/slurm/pairs.tsv"
bash "${BUILD}" "${ANALYSIS}" >"${PAIRS_TSV}"
n_pairs=$(wc -l <"${PAIRS_TSV}" | tr -d ' ')
[[ "${n_pairs}" -gt 0 ]] || {
  echo "ERROR: no pairs generated" >&2
  exit 1
}

# Lighter defaults when the BLAST store is complete; heavier if pairs must BLAST
BLAST_DIR=$(realpath -m "${GENETRIBE_BLAST_DIR:-${WORK}/blast}")
export GENETRIBE_BLAST_DIR="${BLAST_DIR}"
blast_complete=0
if bash "${REPO_ROOT}/scripts/genetribe_check_blasts.sh" "${ANALYSIS}" >/dev/null 2>&1; then
  blast_complete=1
fi
if [[ "${blast_complete}" -eq 1 ]]; then
  GT_CPUS="${GT_CPUS:-16}"
  GT_MEM="${GT_MEM:-64G}"
  GT_TIME="${GT_TIME:-12:00:00}"
  export GENETRIBE_REQUIRE_BLAST="${GENETRIBE_REQUIRE_BLAST:-1}"
else
  GT_CPUS="${GT_CPUS:-48}"
  GT_MEM="${GT_MEM:-128G}"
  GT_TIME="${GT_TIME:-72:00:00}"
  export GENETRIBE_REQUIRE_BLAST="${GENETRIBE_REQUIRE_BLAST:-0}"
  echo "WARN: BLAST store incomplete under ${BLAST_DIR}" >&2
  echo "  Prefer: bash slurm/submit_genetribe_blasts.sh ${GENOME_ARG}" >&2
  echo "  Pair jobs will recompute missing BLASTPs (slow)." >&2
fi

GENETRIBE_SKIP_EXISTING="${GENETRIBE_SKIP_EXISTING:-1}"
: "${GENETRIBE_HOME:=${REPO_ROOT}/genetribe-upstream}"
export GENETRIBE_CONDA_ENV="${GENETRIBE_CONDA_ENV:-}"

if [[ ! -e "${GENETRIBE_HOME}/genetribe" ]]; then
  echo "ERROR: ${GENETRIBE_HOME}/genetribe missing." >&2
  echo "  Conda alone is not enough. From this directory run: make setup" >&2
  echo "  (or: git clone https://github.com/chenym1/genetribe.git genetribe-upstream && cd genetribe-upstream && bash install.sh)" >&2
  exit 1
fi

JOB_SLURM="${REPO_ROOT}/slurm/genetribe_pair_job.slurm"
echo "=== ${n_pairs} GeneTribe pairs (SKIP_EXISTING=${GENETRIBE_SKIP_EXISTING}) ==="
echo "  account=${GT_ACCOUNT:-<default>} partition=${GT_PARTITION:-<default>}"
echo "  per job: CPUS=${GT_CPUS} MEM=${GT_MEM} TIME=${GT_TIME}"
echo "  BLAST_DIR=${BLAST_DIR} complete=${blast_complete} REQUIRE_BLAST=${GENETRIBE_REQUIRE_BLAST}"
echo "  GENETRIBE_HOME=${GENETRIBE_HOME}"
echo "  GENETRIBE_CONDA_ENV=${GENETRIBE_CONDA_ENV:-}"
echo "  GENETRIBE_EXTRA=${GENETRIBE_EXTRA}"
echo "  Slurm stdout/stderr: ${REPO_ROOT}/log/stdout.<pair>.%j.%N"

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
    "${SBATCH_SITE_ARGS[@]}" \
    --job-name="gt-${safe}" \
    --cpus-per-task="${GT_CPUS}" \
    --mem="${GT_MEM}" \
    -t "${GT_TIME}" \
    -N1 \
    --export=ALL,ANALYSIS="${ANALYSIS}",WORK="${WORK}",REPO_ROOT="${REPO_ROOT}",PAIR_L="${a}",PAIR_F="${b}",GENETRIBE_SKIP_EXISTING="${GENETRIBE_SKIP_EXISTING}",GENETRIBE_THREADS="${GT_CPUS}",GENETRIBE_CONDA_ENV="${GENETRIBE_CONDA_ENV:-}",GENETRIBE_HOME="${GENETRIBE_HOME}",GENETRIBE_EXTRA="${GENETRIBE_EXTRA:-}",GENETRIBE_BLAST_DIR="${BLAST_DIR}",GENETRIBE_REQUIRE_BLAST="${GENETRIBE_REQUIRE_BLAST}" \
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
echo "  IDs: ${WORK}/slurm/submitted_job_ids.txt"
echo "Queue: squeue -u \$USER"
echo "When all finish:"
echo "  sbatch slurm/genetribe_finalize.slurm ${GENOME_ARG}"
