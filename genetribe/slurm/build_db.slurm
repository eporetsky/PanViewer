#!/bin/bash
# Build PanViewer *.genetribe.db from subgenome-aware pan TSV.
#
# Usage:
#   sbatch slurm/build_db.slurm wheat
#
# Requires PANVIEWER_ROOT (optional; defaults to parent of genetribe/) and base database/<species>.db.
#
#SBATCH --account=small_grains
#SBATCH --partition=atlas
#SBATCH --job-name=gt-build-db
#SBATCH -N1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH -o "./log/stdout.build_db.%j.%N"
#SBATCH -e "./log/stderr.build_db.%j.%N"

set -euo pipefail

GENOME_ARG="${1:?usage: sbatch slurm/build_db.slurm <species>}"
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then cd "${SLURM_SUBMIT_DIR}" || exit 1; fi
REPO_ROOT="${REPO_ROOT:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}}"
mkdir -p "${REPO_ROOT}/log"


# shellcheck source=/dev/null
source "${REPO_ROOT}/scripts/genetribe_env.sh"
_genetribe_bootstrap || true

bash "${REPO_ROOT}/scripts/genetribe_build_db.sh" "${GENOME_ARG}"
