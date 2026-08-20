#!/usr/bin/env bash
# Build PanViewer annotation DB from subgenome-aware + collinear pan TSV.
#
# Requires:
#   <analysis>/work/genetribe_subgenome_pans.hsh.tsv
#   ${PANVIEWER_ROOT}/database/<species>.db   (base sequences; build_base_index.py)
#   ${PANVIEWER_ROOT}/build_pangene_index.py
#
# Writes:
#   ${PANVIEWER_ROOT}/database/<species>.genetribe.db
#
# Usage:
#   bash scripts/genetribe_build_db.sh wheat
#   PANVIEWER_ROOT=/path/to/panviewer bash scripts/genetribe_build_db.sh oat
set -euo pipefail

GENOME_ARG="${1:?usage: $0 <wheat|oat|barley|…>}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ "${GENOME_ARG}" = /* ]]; then
  ANALYSIS=$(realpath "${GENOME_ARG}")
  SPECIES=$(basename "${ANALYSIS}")
else
  ANALYSIS=$(realpath "${REPO_ROOT}/${GENOME_ARG}")
  SPECIES="${GENOME_ARG}"
fi

# Default: sibling panviewer checkout (ggpanviewer root)
if [[ -z "${PANVIEWER_ROOT:-}" ]]; then
  cand="$(cd "${REPO_ROOT}/.." && pwd)"
  if [[ -f "${cand}/build_pangene_index.py" ]]; then
    PANVIEWER_ROOT="${cand}"
  fi
fi
[[ -n "${PANVIEWER_ROOT:-}" && -d "${PANVIEWER_ROOT}" ]] || {
  echo "ERROR: set PANVIEWER_ROOT to the PanViewer root (dir with build_pangene_index.py)" >&2
  exit 1
}
PANVIEWER_ROOT=$(realpath "${PANVIEWER_ROOT}")

HSH="${ANALYSIS}/work/genetribe_subgenome_pans.hsh.tsv"
BASE="${PANVIEWER_ROOT}/database/${SPECIES}.db"
OUT="${PANVIEWER_ROOT}/database/${SPECIES}.genetribe.db"

[[ -f "${HSH}" ]] || {
  echo "ERROR: missing ${HSH}" >&2
  echo "  Run: python3 scripts/genetribe_cluster_pans.py ${ANALYSIS}" >&2
  exit 1
}
[[ -f "${BASE}" ]] || {
  echo "ERROR: missing base DB ${BASE}" >&2
  echo "  Build sequences first (e.g. python build_base_index.py --force ${SPECIES})" >&2
  exit 1
}
[[ -f "${PANVIEWER_ROOT}/build_pangene_index.py" ]] || {
  echo "ERROR: missing ${PANVIEWER_ROOT}/build_pangene_index.py" >&2
  exit 1
}

echo "=== subgenome pans → ${OUT} ==="
echo "hsh=${HSH}"
echo "base=${BASE}"
ls -lh "${HSH}" "${BASE}"

mkdir -p "${PANVIEWER_ROOT}/pan_sources"
ln -sfn "${HSH}" "${PANVIEWER_ROOT}/pan_sources/${SPECIES}.genetribe.hsh.tsv"

cd "${PANVIEWER_ROOT}"
python3 build_pangene_index.py --force --label genetribe \
  --tsv "${HSH}" --base-db "${BASE}" "${SPECIES}"
ls -lh "${OUT}"
echo "Done: ${OUT}"
