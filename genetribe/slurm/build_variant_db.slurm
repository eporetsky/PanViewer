#!/bin/bash
# Build PanViewer DB for a GeneTribe variant.
#
# Usage:
#   sbatch slurm/build_variant_db.slurm wheat
#
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --job-name=gt-db
#SBATCH -N1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o "./log/stdout.build_db.%j.%N"
#SBATCH -e "./log/stderr.build_db.%j.%N"

set -euo pipefail

GENOME_ARG="${1:?usage: sbatch slurm/build_variant_db.slurm <species>}"
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then cd "${SLURM_SUBMIT_DIR}" || exit 1; fi
REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
mkdir -p "${REPO_ROOT}/log"

bash "${REPO_ROOT}/scripts/genetribe_build_variant_db.sh" "${GENOME_ARG}"
