#!/bin/bash
# RBH connected-component clustering → pan TSV files.
#
# Usage:
#   sbatch slurm/finalize.slurm wheat
#
# Writes:
#   <species>/work/genetribe_pans.{hsh,clust}.tsv
#
# Then build the PanViewer DB:
#   bash scripts/genetribe_build_variant_db.sh wheat
#
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --job-name=gt-finalize
#SBATCH -N1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o "./log/stdout.finalize.%j.%N"
#SBATCH -e "./log/stderr.finalize.%j.%N"

set -euo pipefail

GENOME_ARG="${1:?usage: sbatch slurm/finalize.slurm <species>}"
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then cd "${SLURM_SUBMIT_DIR}" || exit 1; fi
REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
mkdir -p "${REPO_ROOT}/log"

if [[ "${GENOME_ARG}" = /* ]]; then
  ANALYSIS=$(realpath "${GENOME_ARG}")
else
  ANALYSIS=$(realpath "${REPO_ROOT}/${GENOME_ARG}")
fi

WORK=$(realpath -m "${WORK:-${ANALYSIS}/work}")
mkdir -p "${WORK}"

# shellcheck source=/dev/null
source "${REPO_ROOT}/scripts/genetribe_env.sh"
_genetribe_bootstrap || true

PREFIX="${OGI_PREFIX:-OGI}"
python3 "${REPO_ROOT}/scripts/genetribe_rbh_to_pans.py" "${ANALYSIS}" --work "${WORK}" --prefix "${PREFIX}"
