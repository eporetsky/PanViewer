#!/usr/bin/env bash
# Build 05_filtered_pairs.tsv from completed shard DAG outputs (run once before mcl).
# Respects strict_synt from the pan config (same logic as pandagma-pan.sh run_dagchainer).
#
# Usage:
#   pandagma_dagchainer_finalize.sh -c CONFIG -w WORK_DIR
#
set -euo pipefail

CONF=""
WORK_DIR=""

while getopts c:w:h opt; do
  case "$opt" in
    c) CONF=$OPTARG ;;
    w) WORK_DIR=$OPTARG ;;
    h) grep '^#' "$0" | head -12; exit 0 ;;
    *) exit 1 ;;
  esac
done

if [[ -z $CONF || -z $WORK_DIR ]]; then
  echo "ERROR: -c CONFIG -w WORK_DIR required." >&2
  exit 1
fi

CONF=$(realpath "$CONF")
WORK_DIR=$(realpath "$WORK_DIR")
DATA_DIR="${DATA_DIR:-$(dirname "$(dirname "$CONF")")}"

cd "$DATA_DIR" || exit 1
# shellcheck disable=SC1090
source "$CONF"

cd "$WORK_DIR" || exit 1

n_ac=0
shopt -s nullglob
for f in 04_dag/*.aligncoords; do
  [[ $(grep -cve '^#' -e '^$' "$f" 2>/dev/null || echo 0) -gt 0 ]] && n_ac=$((n_ac + 1))
done
shopt -u nullglob

echo "aligncoords files with data: ${n_ac}"
if [[ "${n_ac}" -eq 0 ]]; then
  echo "ERROR: no non-empty 04_dag/*.aligncoords — finish DAG shards first." >&2
  exit 1
fi

if [[ "${strict_synt:-0}" -eq 1 ]]; then
  echo "strict_synt=1: 05_filtered_pairs.tsv from .aligncoords only"
  cat 04_dag/*.aligncoords | awk '$1!~/^#/ {print $2 "\t" $6}' |
    awk 'NF==2' | sort -u > 05_filtered_pairs.tsv
else
  echo "strict_synt=0: 05_filtered_pairs.tsv from matches + .aligncoords"
  cat 04_dag/*_matches.tsv 04_dag/*.aligncoords | awk '$1!~/^#/ {print $2 "\t" $6}' |
    awk 'NF==2' | sort -u > 05_filtered_pairs.tsv
fi

np=$(wc -l < 05_filtered_pairs.tsv | tr -d ' ')
echo "Wrote 05_filtered_pairs.tsv (${np} unique gene pairs)"
