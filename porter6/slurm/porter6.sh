#!/bin/bash
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
#SBATCH --job-name=porter6
#SBATCH --gres=gpu:1
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH -o "./log/stdout.%j.%N"
#SBATCH -e "./log/stderr.%j.%N"

# Submit from the porter6 project root (directory containing porter6/ and slurm/):
#   mkdir -p log
#   sbatch slurm/porter6.sh /path/to/proteins.fasta
#
# Requires: CUDA module, Porter6 installed under porter6/ (see README.md).

set -euo pipefail

_LIB=""
for _c in \
  "${PORTER6_LIB:-}" \
  "${SLURM_SUBMIT_DIR:-}/slurm/porter6.lib.sh" \
  "$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" 2>/dev/null && pwd)/porter6.lib.sh"
do
  [[ -n "${_c:-}" && -f "$_c" ]] && { _LIB="$_c"; break; }
done
unset _c
if [[ -z "$_LIB" ]]; then
  echo "Cannot find porter6.lib.sh. Set PORTER6_LIB or submit from project root." >&2
  exit 1
fi
# shellcheck source=porter6.lib.sh
source "$_LIB"
unset _LIB

FASTA="${1:-}"
if [[ -z "$FASTA" ]]; then
  echo "Usage: sbatch slurm/porter6.sh /path/to/sequences.fasta"
  exit 1
fi

porter6_run_one "$FASTA"
