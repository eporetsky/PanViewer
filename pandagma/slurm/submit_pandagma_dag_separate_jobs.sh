#!/usr/bin/env bash
# Submit N separate Slurm jobs (no job array). Each job: one node, 256GB, 48 CPUs,
# up to 48 parallel run_DAG_chainer.pl processes. Slurm queues excess jobs automatically.
#
# Usage:
#   bash slurm/submit_pandagma_dag_separate_jobs.sh barley
#
# Optional env:
#   DAG_MAX_JOBS=40     never submit more than this (HPC parallel-node quota)
#   DAG_PAIRS_PER_JOB=100  target pairs/job when computing job count
#   DAG_NUM_JOBS=        override job count explicitly
#   DAG_CPUS=48         cpus-per-task per job
#   DAG_NPROC=48        parallel DAG pairs per job (default = DAG_CPUS)
#   DAG_MEM=256GB
#   DAG_TIME=72:00:00
#   DAG_SKIP_EXISTING=1
#
# After all jobs finish:
#   sbatch slurm/pandagma_dagchainer_finalize.slurm barley
#   sbatch slurm/pandagma_pan_resume.slurm barley
#
set -euo pipefail

GENOME_ARG="${1:?usage: $0 <genome e.g. barley>}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

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
mkdir -p "${REPO_ROOT}/log" "${WORK}/slurm"

BUILD="${REPO_ROOT}/scripts/pandagma_dagchainer_build_batches.sh"
chmod +x "${BUILD}" "${REPO_ROOT}/scripts/pandagma_dagchainer_batch_shard.sh" 2>/dev/null || true

DAG_CPUS="${DAG_CPUS:-48}"
DAG_NPROC="${DAG_NPROC:-${DAG_CPUS}}"
DAG_MEM="${DAG_MEM:-256GB}"
DAG_TIME="${DAG_TIME:-72:00:00}"
DAG_SKIP_EXISTING="${DAG_SKIP_EXISTING:-1}"

n_pairs=$(find "${WORK}/04_dag" -maxdepth 1 -name '*_matches.tsv' 2>/dev/null | wc -l | tr -d ' ')
[[ "${n_pairs}" -gt 0 ]] || { echo "ERROR: no 04_dag/*_matches.tsv under ${WORK}" >&2; exit 1; }

DAG_MAX_JOBS="${DAG_MAX_JOBS:-40}"
DAG_PAIRS_PER_JOB="${DAG_PAIRS_PER_JOB:-100}"

if [[ -z "${DAG_NUM_JOBS:-}" ]]; then
  wanted=$(( (n_pairs + DAG_PAIRS_PER_JOB - 1) / DAG_PAIRS_PER_JOB ))
  DAG_NUM_JOBS=$(( wanted < DAG_MAX_JOBS ? wanted : DAG_MAX_JOBS ))
  [[ "${DAG_NUM_JOBS}" -lt 1 ]] && DAG_NUM_JOBS=1
elif [[ "${DAG_NUM_JOBS}" -gt "${DAG_MAX_JOBS}" ]]; then
  echo "WARNING: DAG_NUM_JOBS=${DAG_NUM_JOBS} > DAG_MAX_JOBS=${DAG_MAX_JOBS}; capping at ${DAG_MAX_JOBS}" >&2
  DAG_NUM_JOBS="${DAG_MAX_JOBS}"
fi

pairs_per_job=$(( (n_pairs + DAG_NUM_JOBS - 1) / DAG_NUM_JOBS ))
echo "=== ${n_pairs} match files → ${DAG_NUM_JOBS} separate jobs (max ${DAG_MAX_JOBS}; ~${pairs_per_job} pairs/job) ==="
echo "  per job: CPUS=${DAG_CPUS} NPROC=${DAG_NPROC} MEM=${DAG_MEM} TIME=${DAG_TIME}"
echo "  (dagchainer ≈1 CPU per pair; no job array / no JobArrayTaskLimit)"

bash "${BUILD}" -w "${WORK}" -n "${DAG_NUM_JOBS}"

BATCH_DIR="${WORK}/slurm/dag-batches"
mapfile -t manifests < <(find "${BATCH_DIR}" -maxdepth 1 -name 'batch_*.txt' | sort)
[[ ${#manifests[@]} -gt 0 ]] || { echo "ERROR: no batch_*.txt in ${BATCH_DIR}" >&2; exit 1; }

ENV_FILE="${WORK}/slurm/dag-job-env.sh"
cat > "${ENV_FILE}" <<EOF
REPO_ROOT="${REPO_ROOT}"
CONF="${CONF}"
ANALYSIS="${ANALYSIS}"
WORK="${WORK}"
DAG_NPROC="${DAG_NPROC}"
DAG_SKIP_EXISTING="${DAG_SKIP_EXISTING}"
EOF

JOB_SLURM="${REPO_ROOT}/slurm/pandagma_dagchainer_job.slurm"
submitted=0

for manifest in "${manifests[@]}"; do
  bid=$(basename "${manifest}" .txt)
  jid=$(sbatch --parsable \
    --account=YOUR_ACCOUNT \
    --partition=YOUR_PARTITION \
    --job-name="pdg-dag-${bid}" \
    --cpus-per-task="${DAG_CPUS}" \
    --mem="${DAG_MEM}" \
    -t "${DAG_TIME}" \
    -N1 \
    --export=ALL,WORK="${WORK}",REPO_ROOT="${REPO_ROOT}",BATCH_MANIFEST="${manifest}",BATCH_ID="${bid}" \
    -o "${REPO_ROOT}/log/stdout.${bid}.%j.%N" \
    -e "${REPO_ROOT}/log/stderr.${bid}.%j.%N" \
    "${JOB_SLURM}")
  echo "  ${bid}: job ${jid} ($(wc -l <"${manifest}" | tr -d ' ') pairs)"
  submitted=$((submitted + 1))
done

echo ""
echo "Submitted ${submitted} separate job(s). Queue: squeue -u \$USER"
echo "When all finish:"
echo "  sbatch slurm/pandagma_dagchainer_finalize.slurm ${GENOME_ARG}"
echo "  sbatch slurm/pandagma_pan_resume.slurm ${GENOME_ARG}"
