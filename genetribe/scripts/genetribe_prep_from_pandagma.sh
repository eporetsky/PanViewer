#!/usr/bin/env bash
# Convert Pandagma-style prot/ + bed/ into GeneTribe accessions/<acc>/<acc>.{fa,bed,chrlist}.
#
# Usage:
#   genetribe_prep_from_pandagma.sh <analysis_dir> --chrlist <shared.chrlist>
#   genetribe_prep_from_pandagma.sh <analysis_dir> --prot-dir DIR --bed-dir DIR --chrlist FILE
#
# GeneTribe coreCBS requires **6-column** BED:
#   chrom  start  end  gene_id  score  strand
# Pandagma BEDs are usually **7-column**:
#   chrom  start  end  mRNA_id  score  strand  gene_id
# This script always writes a real 6-col <acc>.bed (never symlink) using gene_id.
#
# Fail-loud: will not invent chromosome groups. Provide --chrlist or place
# accessions/<acc>/<acc>.chrlist yourself.
set -euo pipefail

CHRLIST=""
PROT_DIR=""
BED_DIR=""

usage() {
  cat >&2 <<'EOF'
usage: genetribe_prep_from_pandagma.sh <analysis_dir> [--chrlist FILE]
                                         [--prot-dir DIR] [--bed-dir DIR]
EOF
  exit 1
}

ANALYSIS="${1:?}"
shift || true
ANALYSIS=$(realpath "${ANALYSIS}")

while [[ $# -gt 0 ]]; do
  case "$1" in
    --chrlist) CHRLIST=$(realpath "$2"); shift 2 ;;
    --prot-dir) PROT_DIR=$(realpath "$2"); shift 2 ;;
    --bed-dir) BED_DIR=$(realpath "$2"); shift 2 ;;
    -h|--help) usage ;;
    *) echo "ERROR: unknown arg $1" >&2; usage ;;
  esac
done

PROT_DIR="${PROT_DIR:-${ANALYSIS}/prot}"
BED_DIR="${BED_DIR:-${ANALYSIS}/bed}"

[[ -d "${PROT_DIR}" ]] || { echo "ERROR: missing prot dir ${PROT_DIR}" >&2; exit 1; }
[[ -d "${BED_DIR}" ]] || { echo "ERROR: missing bed dir ${BED_DIR}" >&2; exit 1; }

# Write GeneTribe 6-col BED from 6- or 7-col input (stdout).
to_genetribe_bed6() {
  local src="$1"
  local opener=(cat)
  [[ "${src}" == *.gz ]] && opener=(gzip -dc)
  "${opener[@]}" "${src}" | gawk -vOFS='\t' '
    /^#/ || NF == 0 { next }
    NF == 6 {
      print $1, $2, $3, $4, $5, $6
      next
    }
    NF >= 7 {
      # Pandagma: gene id in column 7
      print $1, $2, $3, $7, $5, $6
      next
    }
    {
      printf("ERROR: %s line %d has %d columns (need 6 or 7)\n", FILENAME, NR, NF) > "/dev/stderr"
      exit 1
    }
  '
}

n=0
for fa in "${PROT_DIR}"/*.fa "${PROT_DIR}"/*.faa "${PROT_DIR}"/*.fasta; do
  [[ -e "${fa}" ]] || continue
  base=$(basename "${fa}")
  acc="${base%.*}"
  bed=""
  for cand in "${BED_DIR}/${acc}.bed" "${BED_DIR}/${acc}.bed6" "${BED_DIR}/${acc}.bed.gz"; do
    if [[ -f "${cand}" ]]; then
      bed="${cand}"
      break
    fi
  done
  [[ -n "${bed}" ]] || {
    echo "ERROR: no BED for accession '${acc}' under ${BED_DIR}" >&2
    exit 1
  }

  dest="${ANALYSIS}/accessions/${acc}"
  mkdir -p "${dest}"
  ln -sfn "${fa}" "${dest}/${acc}.fa"

  # Always materialize 6-col bed (break any old symlink to 7-col Pandagma bed)
  rm -f "${dest}/${acc}.bed"
  to_genetribe_bed6 "${bed}" >"${dest}/${acc}.bed"
  ncols=$(gawk 'NF && $1 !~ /^#/ { print NF; exit }' "${dest}/${acc}.bed")
  [[ "${ncols}" == "6" ]] || {
    echo "ERROR: ${dest}/${acc}.bed has ${ncols:-0} columns after conversion (need 6)" >&2
    exit 1
  }

  if [[ -n "${CHRLIST}" ]]; then
    cp -f "${CHRLIST}" "${dest}/${acc}.chrlist"
  elif [[ ! -f "${dest}/${acc}.chrlist" ]]; then
    echo "ERROR: missing ${dest}/${acc}.chrlist (pass --chrlist or create it)" >&2
    exit 1
  fi
  n=$((n + 1))
done

[[ "${n}" -gt 0 ]] || {
  echo "ERROR: no protein FASTA files under ${PROT_DIR}" >&2
  exit 1
}
echo "Prepared ${n} accession(s) under ${ANALYSIS}/accessions/ (6-col GeneTribe BEDs)"
