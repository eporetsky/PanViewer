#!/usr/bin/env bash
# Run DAGchainer for all pairs listed in one batch manifest (one separate Slurm job).
#
# Usage:
#   pandagma_dagchainer_batch_shard.sh -c CONFIG -d DATA_DIR -w WORK_DIR -m MANIFEST
#
set -euo pipefail

CONF=""
DATA_DIR=""
WORK_DIR=""
MANIFEST=""

while getopts c:d:w:m:h opt; do
  case "$opt" in
    c) CONF=$OPTARG ;;
    d) DATA_DIR=$OPTARG ;;
    w) WORK_DIR=$OPTARG ;;
    m) MANIFEST=$OPTARG ;;
    h) grep '^#' "$0" | head -12; exit 0 ;;
    *) exit 1 ;;
  esac
done

if [[ -z $CONF || -z $DATA_DIR || -z $WORK_DIR || -z $MANIFEST ]]; then
  echo "ERROR: -c CONFIG -d DATA_DIR -w WORK_DIR -m MANIFEST required." >&2
  exit 1
fi

CONF=$(realpath "$CONF")
DATA_DIR=$(realpath "$DATA_DIR")
WORK_DIR=$(realpath "$WORK_DIR")
MANIFEST=$(realpath "$MANIFEST")

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=pandagma_dagchainer_lib.sh
source "${REPO_ROOT}/scripts/pandagma_dagchainer_lib.sh"
pandagma_dagchainer_load_config || exit 1

mapfile -t matches < <(grep -v '^[[:space:]]*$' "${MANIFEST}")
if [[ ${#matches[@]} -eq 0 ]]; then
  echo "ERROR: empty manifest ${MANIFEST}" >&2
  exit 1
fi

echo "DAG batch $(basename "${MANIFEST}") pairs=${#matches[@]} DAG_NPROC=${DAG_NPROC}"
pandagma_dagchainer_run_pairs "${matches[@]}"
