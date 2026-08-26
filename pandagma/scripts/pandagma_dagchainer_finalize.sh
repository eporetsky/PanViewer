#!/usr/bin/env bash
# Build 05_filtered_pairs.tsv from completed shard DAG outputs (run once before mcl).
# Respects strict_synt from the pan config (same logic as pandagma-pan.sh run_dagchainer).
#
# Large shard counts: per-file pair extraction runs in parallel (GNU parallel when
# available), each shard is sort -u'd locally, then chunks are merged with sort -u -m
# (GNU sort --parallel when available).
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
    h) grep '^#' "$0" | head -16; exit 0 ;;
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

n_cpus() {
  if command -v nproc &>/dev/null; then
    nproc
  elif [[ -r /proc/cpuinfo ]]; then
    grep -c ^processor /proc/cpuinfo
  else
    sysctl -n hw.ncpu 2>/dev/null || echo 4
  fi
}

# Prefer GNU sort (Linux or brew gsort on macOS) for --parallel / -m merge.
sort_cmd() {
  if command -v gsort &>/dev/null && gsort --version 2>/dev/null | grep -q GNU; then
    echo gsort
  elif sort --version 2>/dev/null | grep -q GNU; then
    echo sort
  else
    echo sort
  fi
}

sort_is_gnu() {
  local s
  s=$(sort_cmd)
  "$s" --version 2>/dev/null | grep -q GNU
}

sort_parallel_flags() {
  if sort_is_gnu; then
    echo --parallel="$(n_cpus)" -S 50%
  fi
}

# At least one non-comment, non-empty line.
file_has_data() {
  awk 'NF && $1 !~ /^#/ { exit 0 } END { exit 1 }' "$1"
}

extract_pairs_chunk() {
  local infile=$1 outchunk=$2 sort_bin=$3 work_dir=$4
  LC_ALL=C awk '$1 !~ /^#/ && NF >= 6 { print $2 "\t" $6 }' "$infile" |
    "$sort_bin" -u -T "$work_dir" -o "$outchunk"
}

cd "$DATA_DIR" || exit 1
# shellcheck disable=SC1090
source "$CONF"

cd "$WORK_DIR" || exit 1

NCPUS=$(n_cpus)
SORT_BIN=$(sort_cmd)
# shellcheck disable=SC2206
SORT_PAR=($(sort_parallel_flags))

shopt -s nullglob
aligncoords=(04_dag/*.aligncoords)
if ((${#aligncoords[@]} == 0)); then
  echo "ERROR: no 04_dag/*.aligncoords — finish DAG shards first." >&2
  exit 1
fi

n_ac=0
if command -v parallel &>/dev/null; then
  n_ac=$(
    parallel -j "$NCPUS" --halt soon,fail=1 '
      awk "NF && \$1 !~ /^#/" {} >/dev/null && echo 1
    ' ::: "${aligncoords[@]}" | wc -l | tr -d ' '
  )
else
  for f in "${aligncoords[@]}"; do
    file_has_data "$f" && n_ac=$((n_ac + 1))
  done
fi

echo "aligncoords files with data: ${n_ac}"
if [[ "${n_ac}" -eq 0 ]]; then
  echo "ERROR: no non-empty 04_dag/*.aligncoords — finish DAG shards first." >&2
  exit 1
fi

if [[ "${strict_synt:-0}" -eq 1 ]]; then
  echo "strict_synt=1: 05_filtered_pairs.tsv from .aligncoords only"
  inputs=("${aligncoords[@]}")
else
  echo "strict_synt=0: 05_filtered_pairs.tsv from matches + .aligncoords"
  inputs=(04_dag/*_matches.tsv "${aligncoords[@]}")
fi

CHUNK_DIR=$(mktemp -d "${WORK_DIR}/05_pair_chunks.XXXXXX")
cleanup() { rm -rf "$CHUNK_DIR"; }
trap cleanup EXIT

echo "Extracting pairs (${#inputs[@]} files, ${NCPUS} jobs, sort=${SORT_BIN})"

if command -v parallel &>/dev/null && ((${#inputs[@]} > 1)); then
  export SORT_BIN WORK_DIR CHUNK_DIR
  parallel -j "$NCPUS" --halt soon,fail=1 '
    LC_ALL=C awk '"'"'$1 !~ /^#/ && NF >= 6 { print $2 "\t" $6 }'"'"' {} |
      '"$SORT_BIN"' -u -T "'"$WORK_DIR"'" -o "'"$CHUNK_DIR"'/$(basename {}).pairs"
  ' ::: "${inputs[@]}"
else
  for f in "${inputs[@]}"; do
    extract_pairs_chunk "$f" "${CHUNK_DIR}/$(basename "$f").pairs" "$SORT_BIN" "$WORK_DIR"
  done
fi

chunk_files=("$CHUNK_DIR"/*.pairs)
if ((${#chunk_files[@]} == 0)); then
  echo "ERROR: no pair chunks produced." >&2
  exit 1
fi

echo "Merging ${#chunk_files[@]} sorted chunks -> 05_filtered_pairs.tsv"
# shellcheck disable=SC2086
LC_ALL=C "$SORT_BIN" -u -m ${SORT_PAR[@]+"${SORT_PAR[@]}"} -T . "${chunk_files[@]}" -o 05_filtered_pairs.tsv

np=$(wc -l < 05_filtered_pairs.tsv | tr -d ' ')
echo "Wrote 05_filtered_pairs.tsv (${np} unique gene pairs)"
