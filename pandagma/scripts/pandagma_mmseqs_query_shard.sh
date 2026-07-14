#!/usr/bin/env bash
# Run one "query" slice of Pandagma's pan mmseqs step: pairwise mmseqs easy-cluster for
# every pair (QUERY_INDEX, j) with j > QUERY_INDEX in cds_files order (upper triangle).
#
# Same outputs as pandagma-upstream/bin/pandagma-pan.sh run_mmseqs():
#   WORK_DIR/03_mmseqs/${qry}.x.${sbj}_cluster.tsv
#
# Prerequisites: ingest completed (WORK_DIR/02_fasta_nuc/*.${fa} present).
#
# Usage:
#   pandagma_mmseqs_query_shard.sh -c CONFIG -d DATA_DIR -w WORK_DIR -i QUERY_INDEX [-n THREADS] [-s] [-f]
#
#   -n THREADS   mmseqs --threads (default: 4)
#   -s           Skip pair if *_cluster.tsv already exists
#   -f           After mmseqs, run member-expansion filter on this shard's cluster files only
#                (03_mmseqs/<query>.x.*_cluster.tsv -> 04_dag/<stem>_matches.tsv)

set -euo pipefail

usage() {
  grep '^#' "$0" | grep -v '^#!/' | sed 's/^# \{0,1\}//' | head -26
  exit "${1:-0}"
}

CONF=""
DATA_DIR=""
WORK_DIR=""
QUERY_INDEX=""
THREADS=4
SKIP_EXISTING=0
RUN_FILTER=0

while getopts c:d:w:i:n:sfh opt; do
  case "$opt" in
    c) CONF=$OPTARG ;;
    d) DATA_DIR=$OPTARG ;;
    w) WORK_DIR=$OPTARG ;;
    i) QUERY_INDEX=$OPTARG ;;
    n) THREADS=$OPTARG ;;
    s) SKIP_EXISTING=1 ;;
    f) RUN_FILTER=1 ;;
    h) usage 0 ;;
    *) usage 1 ;;
  esac
done

if [[ -z $CONF || -z $DATA_DIR || -z $WORK_DIR || -z $QUERY_INDEX ]]; then
  echo "ERROR: -c CONFIG -d DATA_DIR -w WORK_DIR -i QUERY_INDEX are required." >&2
  usage 1
fi

if ! [[ "$QUERY_INDEX" =~ ^[0-9]+$ ]]; then
  echo "ERROR: -i must be a non-negative integer (got: $QUERY_INDEX)" >&2
  exit 1
fi

CONF=$(realpath "$CONF")
DATA_DIR=$(realpath "$DATA_DIR")
WORK_DIR=$(realpath "$WORK_DIR")

if [[ ! -d $WORK_DIR/02_fasta_nuc ]]; then
  echo "ERROR: $WORK_DIR/02_fasta_nuc not found. Run pandagma pan ... -s ingest first." >&2
  exit 1
fi

declare -a cds_files
cd "$DATA_DIR" || exit 1
# shellcheck disable=SC1090
source "$CONF"
mapfile -t cds_files < <(realpath --canonicalize-existing "${cds_files[@]}")

n=${#cds_files[@]}
if (( n < 2 )); then
  echo "ERROR: cds_files must list at least two FASTA paths in the config." >&2
  exit 1
fi

if (( QUERY_INDEX >= n )); then
  echo "ERROR: QUERY_INDEX ($QUERY_INDEX) must be < number of genomes ($n)." >&2
  exit 1
fi

fasta_file=$(basename "${cds_files[0]}" .gz)
fa="${fasta_file##*.}"

file1_num=$QUERY_INDEX
if (( file1_num + 1 >= n )); then
  echo "Query index $file1_num is the last genome in cds_files order: no pairs (i,j) with j>i. Nothing to do."
  exit 0
fi

mkdir -p "$WORK_DIR/03_mmseqs" "$WORK_DIR/03_mmseqs_tmp"
cd "$WORK_DIR" || exit 1

qry_base=$(basename "${cds_files[file1_num]%.*}" ."$fa")

for (( file2_num = file1_num + 1; file2_num < n; file2_num++ )); do
  sbj_base=$(basename "${cds_files[file2_num]%.*}" ."$fa")
  outbase="${qry_base}.x.${sbj_base}"
  cluster_tsv="03_mmseqs/${outbase}_cluster.tsv"

  if [[ $SKIP_EXISTING -eq 1 && -f $cluster_tsv ]]; then
    echo "  Skip existing: $cluster_tsv"
    continue
  fi

  echo "  Running mmseqs on comparison: ${qry_base}.x.${sbj_base}"
  MMTEMP=$(mktemp -d -p 03_mmseqs_tmp)
  # Pandagma sets MMSEQS_NUM_THREADS in main_pan_fam; mmseqs respects it if CLI threads omitted.
  export MMSEQS_NUM_THREADS="$THREADS"
  {
    cat "02_fasta_nuc/${qry_base}.${fa}" "02_fasta_nuc/${sbj_base}.${fa}"
  } | mmseqs easy-cluster stdin "03_mmseqs/${outbase}" "$MMTEMP" \
      --min-seq-id "$clust_iden" -c "$clust_cov" --cov-mode 0 --cluster-reassign 1>/dev/null

  rm -f "03_mmseqs/${outbase}_rep_seq.fasta" "03_mmseqs/${outbase}_all_seqs.fasta"
done

if [[ $RUN_FILTER -eq 1 ]]; then
  EXPECTED_INC=""
  shopt -s nullglob
  for cand in "${DATA_DIR}/config/"*expected_chr_matches*.inc; do
    EXPECTED_INC="$cand"
    break
  done
  shopt -u nullglob
  if [[ -z $EXPECTED_INC ]]; then
    echo "ERROR: -f requires an include under ${DATA_DIR}/config/*expected_chr_matches*.inc" >&2
    exit 1
  fi
  REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
  FILTER_PY="${REPO_ROOT}/scripts/pandagma_mmseqs_filter.py"
  if [[ ! -f "$FILTER_PY" ]]; then
    echo "ERROR: filter script missing: $FILTER_PY" >&2
    exit 1
  fi
  echo "Shard filter (--expected-inc from config): expanding 03_mmseqs/${qry_base}.x.*_cluster.tsv"
  shopt -s nullglob
  for cluster in "03_mmseqs/${qry_base}".x.*_cluster.tsv; do
    echo "  $cluster"
    python3 "$FILTER_PY" --cluster-tsv "$(pwd)/$cluster" --expected-inc "$(realpath "$EXPECTED_INC")"
  done
  shopt -u nullglob
fi

echo "Done mmseqs shard for query index $QUERY_INDEX ($qry_base)."
