#!/usr/bin/env bash
# Run on the *login node* (or another shell), not via sbatch.
# Queues one GPU job per FASTA; jobs can wait in the queue as long as needed.
#
#   cd /path/to/project    # contains porter6/ and slurm/
#   bash slurm/porter6.batch.sh /path/to/dir_of_fastas
#
# Skips FASTAs that already have results/<stem>/ (directory exists).
# Each job sets PORTER6_ISOLATE=1 so it uses a private copy of porter6/ (see porter6.lib.sh).

set -euo pipefail

_SLURM_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${_SLURM_DIR}/.." && pwd)"
# shellcheck source=porter6.lib.sh
source "${_SLURM_DIR}/porter6.lib.sh"

FASTA_DIR="${1:-}"
if [[ -z "$FASTA_DIR" || ! -d "$FASTA_DIR" ]]; then
  echo "Usage: bash slurm/porter6.batch.sh /path/to/dir_with_fastas"
  exit 1
fi

cd "$ROOT"

if [[ "$FASTA_DIR" != /* ]]; then
  FASTA_DIR="${ROOT}/${FASTA_DIR}"
fi

mapfile -t _FASTAS < <(find "$FASTA_DIR" -maxdepth 1 -type f \( -iname '*.fasta' -o -iname '*.fa' -o -iname '*.faa' \) | LC_ALL=C sort)

if [[ ${#_FASTAS[@]} -eq 0 ]]; then
  echo "No FASTA files found in: $FASTA_DIR"
  exit 1
fi

echo "Submitting up to ${#_FASTAS[@]} job(s) from $FASTA_DIR"

for FASTA in "${_FASTAS[@]}"; do
  STEM=$(porter6_stem_from_fasta "$FASTA")
  RESULTS_DIR="${ROOT}/results/${STEM}"
  if [[ -d "$RESULTS_DIR" ]]; then
    echo "Skipping $(basename "$FASTA") (already exists: $RESULTS_DIR)"
    continue
  fi
  echo "Queueing $(basename "$FASTA") -> ${RESULTS_DIR}"
  # Pass absolute path so compute jobs find the lib on shared disk (spool copy breaks BASH_SOURCE).
  sbatch --export=ALL,PORTER6_ISOLATE=1,PORTER6_LIB="${_SLURM_DIR}/porter6.lib.sh" "${_SLURM_DIR}/porter6.sh" "$FASTA"
done

echo "Done submitting."
