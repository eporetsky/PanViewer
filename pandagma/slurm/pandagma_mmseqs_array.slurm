#!/bin/bash
# Submit one mmseqs job per genome in cds_files (same order as config pan conf).
# Needs ingest done first (02_fasta_nuc/ under work dir).
# After each shard's mmseqs pairs finish, that job runs member-expansion filter on its cluster
# files only (writes 04_dag/<stem>_matches.tsv). No separate bulk filter Slurm step.
# Usage: sbatch slurm/pandagma_mmseqs_array.slurm wheat
# Optional env: CONF (absolute path), WORK, PANDAGMA_CONDA_ENV (conda env with mmseqs).

#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --job-name=pdg-mmseqs-launch
#SBATCH -N1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH -t 00:30:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

set -euo pipefail

GENOME_ARG="${1:?usage: sbatch slurm/pandagma_mmseqs_array.slurm <genome e.g. wheat>}"
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then cd "${SLURM_SUBMIT_DIR}" || exit 1; fi
REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"

if [[ "${GENOME_ARG}" = /* ]]; then
  ANALYSIS=$(realpath "${GENOME_ARG}")
else
  ANALYSIS=$(realpath "${REPO_ROOT}/${GENOME_ARG}")
fi

if [[ ! -d "${ANALYSIS}/cds" ]] || [[ ! -d "${ANALYSIS}/prot" ]]; then
  echo "ERROR: need ${ANALYSIS}/cds and ${ANALYSIS}/prot (preprocess inputs first)." >&2
  exit 1
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
mkdir -p "${ANALYSIS}/work/slurm" "${WORK}" "${REPO_ROOT}/log"

LIST="${REPO_ROOT}/scripts/pandagma_list_cds_basenames.sh"
[[ -x "${LIST}" ]] || chmod +x "${LIST}" 2>/dev/null || true

n=0
while IFS=$'\t' read -r idx gid; do
  [[ -z "${idx:-}" ]] && continue
  gid="${gid//$'\r'/}"
  safe="${gid//[^A-Za-z0-9._-]/_}"
  out_sl="${ANALYSIS}/work/slurm/pdg-mmseqs-${safe}.slurm"
  cat > "${out_sl}" <<EOF
#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --job-name=pdg-mmseqs
#SBATCH -N1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=256GB
#SBATCH -t 72:00:00
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"
set -euo pipefail
REPO_ROOT="${REPO_ROOT}"
export PATH="\${REPO_ROOT}/pandagma-upstream/bin:\$PATH"
[[ -n "\${PANDAGMA_CONDA_ENV:-}" ]] && export PATH="\${PANDAGMA_CONDA_ENV}/bin:\$PATH"
export SLURM_CPUS_PER_TASK="\${SLURM_CPUS_PER_TASK:-48}"
bash "\${REPO_ROOT}/scripts/pandagma_mmseqs_query_shard.sh" \\
  -c "${CONF}" -d "${ANALYSIS}" -w "${WORK}" -i "${idx}" -n "\${SLURM_CPUS_PER_TASK}" -s -f
EOF
  chmod +x "${out_sl}"
  sbatch "${out_sl}"
  n=$((n + 1))
done < <(bash "${LIST}" -c "${CONF}" -d "${ANALYSIS}")

[[ "${n}" -gt 0 ]] || { echo "ERROR: no genomes from ${LIST}" >&2; exit 1; }
echo "Submitted ${n} job(s); each shard runs mmseqs then filter for its pairwise outputs. conf=${CONF} work=${WORK}"
