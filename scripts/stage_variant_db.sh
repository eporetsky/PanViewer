#!/usr/bin/env bash
# Stage a pan-gene variant under input/ and build database/<variant_id>.db
#
# Usage:
#   bash scripts/stage_variant_db.sh wheat.genetribe \
#     --pan-tsv genetribe/wheat/work/genetribe_pans.hsh.tsv \
#     --bed-dir genetribe/wheat/bed \
#     --prot-dir genetribe/wheat/prot
#
# Optional:
#   --cds-dir PATH       CDS FASTA directory (symlinked if present)
#   --porter6-dir PATH   Porter6 CSV directory (symlinked if present)
#   --dry-run            Print plan only
#   --no-force           Skip --force on build_index.py
#
# The staging id is the input/ folder name and database/<id>.db stem (e.g. wheat.pandagma).
# Map it to a UI variant in database/config.json — see database/config.json.example.
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VARIANT_ID=""
PAN_TSV=""
BED_DIR=""
PROT_DIR=""
CDS_DIR=""
PORTER6_DIR=""
DRY_RUN=0
FORCE=1

usage() {
  cat >&2 <<EOF
usage: $0 <variant_id> --pan-tsv PATH --bed-dir PATH --prot-dir PATH [options]

  variant_id   DB stem, e.g. wheat.pandagma, wheat.genetribe, oat.panoat, barley.panbarlex
  --pan-tsv      Pandagma/GeneTribe/OrthoFinder pan membership (*.hsh.tsv or *.clust.tsv)
  --bed-dir      BED directory (one file per accession)
  --prot-dir     Protein FASTA directory

Options: --cds-dir, --porter6-dir, --dry-run, --no-force
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pan-tsv) PAN_TSV="$2"; shift 2 ;;
    --bed-dir) BED_DIR="$2"; shift 2 ;;
    --prot-dir) PROT_DIR="$2"; shift 2 ;;
    --cds-dir) CDS_DIR="$2"; shift 2 ;;
    --porter6-dir) PORTER6_DIR="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1; shift ;;
    --no-force) FORCE=0; shift ;;
    -h|--help) usage ;;
    -*)
      echo "ERROR: unknown option: $1" >&2
      usage
      ;;
    *)
      if [[ -z "${VARIANT_ID}" ]]; then
        VARIANT_ID="$1"
      else
        echo "ERROR: unexpected argument: $1" >&2
        usage
      fi
      shift
      ;;
  esac
done

[[ -n "${VARIANT_ID}" && -n "${PAN_TSV}" && -n "${BED_DIR}" && -n "${PROT_DIR}" ]] || usage

VARIANT_ID="$(echo "${VARIANT_ID}" | tr '[:upper:]' '[:lower:]')"
PAN_TSV="$(realpath "${PAN_TSV}")"
BED_DIR="$(realpath "${BED_DIR}")"
PROT_DIR="$(realpath "${PROT_DIR}")"
[[ -n "${CDS_DIR}" ]] && CDS_DIR="$(realpath "${CDS_DIR}")"
[[ -n "${PORTER6_DIR}" ]] && PORTER6_DIR="$(realpath "${PORTER6_DIR}")"

INPUT_DIR="${REPO_ROOT}/input/${VARIANT_ID}"
DB_PATH="${REPO_ROOT}/database/${VARIANT_ID}.db"

for need in "${PAN_TSV}" "${BED_DIR}" "${PROT_DIR}"; do
  [[ -e "${need}" ]] || { echo "ERROR: not found: ${need}" >&2; exit 1; }
done

echo "=== stage_variant_db ==="
echo "  variant=${VARIANT_ID}"
echo "  input=${INPUT_DIR}"
echo "  db=${DB_PATH}"

stage() {
  mkdir -p "${INPUT_DIR}" "${REPO_ROOT}/database"
  local pan_name
  pan_name="$(basename "${PAN_TSV}")"
  ln -sfn "${PAN_TSV}" "${INPUT_DIR}/${pan_name}"
  ln -sfn "${BED_DIR}" "${INPUT_DIR}/bed"
  ln -sfn "${PROT_DIR}" "${INPUT_DIR}/prot"
  [[ -n "${CDS_DIR}" && -e "${CDS_DIR}" ]] && ln -sfn "${CDS_DIR}" "${INPUT_DIR}/cds"
  [[ -n "${PORTER6_DIR}" && -e "${PORTER6_DIR}" ]] && ln -sfn "${PORTER6_DIR}" "${INPUT_DIR}/porter6"
}

if [[ "${DRY_RUN}" -eq 1 ]]; then
  echo "dry-run: would stage ${INPUT_DIR} and run build_index.py ${VARIANT_ID}"
  exit 0
fi

stage
cd "${REPO_ROOT}"
args=()
[[ "${FORCE}" -eq 1 ]] && args+=(--force)
args+=("${VARIANT_ID}")
python build_index.py "${args[@]}"
ls -lh "${DB_PATH}"
echo "Done: ${DB_PATH}"
