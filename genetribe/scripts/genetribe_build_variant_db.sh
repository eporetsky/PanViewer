#!/usr/bin/env bash
# Build PanViewer SQLite for a GeneTribe variant (wrapper around stage_variant_db.sh).
#
# Usage:
#   bash scripts/genetribe_build_variant_db.sh wheat
#
# Env overrides:
#   VARIANT_ID=wheat.genetribe
#   HSH=/path/to/genetribe_pans.hsh.tsv
#
set -euo pipefail

GENOME_ARG="${1:?usage: $0 <analysis e.g. wheat>}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PANVIEWER_ROOT="$(cd "${REPO_ROOT}/.." && pwd)"

if [[ "${GENOME_ARG}" = /* ]]; then
  ANALYSIS=$(realpath "${GENOME_ARG}")
  SPECIES=$(basename "${ANALYSIS}")
else
  ANALYSIS=$(realpath "${REPO_ROOT}/${GENOME_ARG}")
  SPECIES="${GENOME_ARG}"
fi

VARIANT_ID="${VARIANT_ID:-${SPECIES}.genetribe}"
HSH="${HSH:-${ANALYSIS}/work/genetribe_pans.hsh.tsv}"

[[ -f "${HSH}" ]] || {
  echo "ERROR: missing ${HSH}. Run: sbatch slurm/finalize.slurm ${SPECIES}" >&2
  exit 1
}

ARGS=(
  "${VARIANT_ID}"
  --pan-tsv "${HSH}"
  --bed-dir "${ANALYSIS}/bed"
  --prot-dir "${ANALYSIS}/prot"
)
[[ -d "${ANALYSIS}/cds" || -L "${ANALYSIS}/cds" ]] && ARGS+=(--cds-dir "${ANALYSIS}/cds")
[[ -d "${ANALYSIS}/porter6" || -L "${ANALYSIS}/porter6" ]] && ARGS+=(--porter6-dir "${ANALYSIS}/porter6")

bash "${PANVIEWER_ROOT}/scripts/stage_variant_db.sh" "${ARGS[@]}"
