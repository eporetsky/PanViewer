#!/usr/bin/env bash
# Print CDS FASTA basenames in the same order as pandagma uses for cds_files
# (after stripping extensions like pandagma-pan.sh run_mmseqs). For sanity checks
# against SLURM_ARRAY_TASK_ID when using pandagma_mmseqs_query_shard.sh.
#
# Usage: pandagma_list_cds_basenames.sh -c CONFIG -d DATA_DIR

set -euo pipefail

CONF=""
DATA_DIR=""

while getopts c:d:h opt; do
  case "$opt" in
    c) CONF=$OPTARG ;;
    d) DATA_DIR=$OPTARG ;;
    h) head -20 "$0" | tail -19; exit 0 ;;
    *) exit 1 ;;
  esac
done

if [[ -z $CONF || -z $DATA_DIR ]]; then
  echo "Usage: $0 -c CONFIG -d DATA_DIR" >&2
  exit 1
fi

CONF=$(realpath "$CONF")
cd "$(realpath "$DATA_DIR")" || exit 1
declare -a cds_files
# shellcheck disable=SC1090
source "$CONF"
mapfile -t cds_files < <(realpath --canonicalize-existing "${cds_files[@]}")

fasta_file=$(basename "${cds_files[0]}" .gz)
fa="${fasta_file##*.}"

i=0
for _ in "${cds_files[@]}"; do
  base=$(basename "${cds_files[i]%.*}" ."$fa")
  printf '%d\t%s\n' "$i" "$base"
  i=$((i + 1))
done
