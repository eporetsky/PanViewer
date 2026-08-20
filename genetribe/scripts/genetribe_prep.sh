#!/usr/bin/env bash
# Prepare GeneTribe accession inputs from protein FASTA + BED or GFF.
#
# Layout written:
#   <analysis>/accessions/<acc>/<acc>.{fa,bed,chrlist}
#
# Usage:
#   # Already have BED (6-col GeneTribe or 7-col Pandagma):
#   bash scripts/genetribe_prep.sh wheat --chrlist examples/wheat.chrlist \
#     --prot-dir wheat/prot --bed-dir wheat/bed
#
#   # Protein FASTA + GFF/GFF3 (gene features → 6-col BED):
#   bash scripts/genetribe_prep.sh wheat --chrlist examples/wheat.chrlist \
#     --prot-dir wheat/prot --gff-dir wheat/gff
#
# Defaults: --prot-dir <analysis>/prot, --bed-dir <analysis>/bed
#           (or --gff-dir <analysis>/gff if that directory exists and bed does not)
#
# Gene IDs in FASTA headers must match BED column 4 (gene-level IDs only).
# Fail-loud: will not invent chromosome groups — pass --chrlist.
set -euo pipefail

CHRLIST=""
PROT_DIR=""
BED_DIR=""
GFF_DIR=""

usage() {
  cat >&2 <<'EOF'
usage: genetribe_prep.sh <analysis_dir> --chrlist FILE
                         [--prot-dir DIR] [--bed-dir DIR | --gff-dir DIR]
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
    --gff-dir) GFF_DIR=$(realpath "$2"); shift 2 ;;
    -h|--help) usage ;;
    *) echo "ERROR: unknown arg $1" >&2; usage ;;
  esac
done

PROT_DIR="${PROT_DIR:-${ANALYSIS}/prot}"
[[ -d "${PROT_DIR}" ]] || { echo "ERROR: missing prot dir ${PROT_DIR}" >&2; exit 1; }
[[ -n "${CHRLIST}" && -f "${CHRLIST}" ]] || {
  echo "ERROR: --chrlist FILE is required (GeneTribe chromosome-group patterns)" >&2
  exit 1
}

if [[ -n "${BED_DIR}" && -n "${GFF_DIR}" ]]; then
  echo "ERROR: pass only one of --bed-dir or --gff-dir" >&2
  exit 1
fi
if [[ -z "${BED_DIR}" && -z "${GFF_DIR}" ]]; then
  if [[ -d "${ANALYSIS}/bed" ]]; then
    BED_DIR=$(realpath "${ANALYSIS}/bed")
  elif [[ -d "${ANALYSIS}/gff" ]]; then
    GFF_DIR=$(realpath "${ANALYSIS}/gff")
  else
    echo "ERROR: need --bed-dir or --gff-dir (or <analysis>/bed or <analysis>/gff)" >&2
    exit 1
  fi
fi

AWK=(awk)
command -v gawk >/dev/null 2>&1 && AWK=(gawk)

# Write GeneTribe 6-col BED from 6- or 7-col input (stdout).
to_genetribe_bed6() {
  local src="$1"
  local opener=(cat)
  [[ "${src}" == *.gz ]] && opener=(gzip -dc)
  "${opener[@]}" "${src}" | "${AWK[@]}" -vOFS='\t' '
    /^#/ || NF == 0 { next }
    NF == 6 {
      print $1, $2, $3, $4, $5, $6
      next
    }
    NF >= 7 {
      print $1, $2, $3, $7, $5, $6
      next
    }
    {
      printf("ERROR: %s line %d has %d columns (need 6 or 7)\n", FILENAME, NR, NF) > "/dev/stderr"
      exit 1
    }
  '
}

# gene features from GFF/GFF3 → 6-col BED (ID / Name / gene_id attributes).
gff_to_bed6() {
  local src="$1"
  local opener=(cat)
  [[ "${src}" == *.gz ]] && opener=(gzip -dc)
  "${opener[@]}" "${src}" | "${AWK[@]}" -vOFS='\t' '
    BEGIN { FS = "\t" }
    /^#/ || NF < 9 { next }
    $3 != "gene" { next }
    {
      chrom = $1; start = $4 - 1; end = $5; strand = $7
      if (strand != "+" && strand != "-") strand = "."
      attrs = $9
      gid = ""
      n = split(attrs, parts, ";")
      for (i = 1; i <= n; i++) {
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", parts[i])
        if (parts[i] ~ /^ID=/) { gid = substr(parts[i], 4); break }
      }
      if (gid == "") {
        for (i = 1; i <= n; i++) {
          if (parts[i] ~ /^gene_id=/) { gid = substr(parts[i], 9); break }
        }
      }
      if (gid == "") {
        for (i = 1; i <= n; i++) {
          if (parts[i] ~ /^Name=/) { gid = substr(parts[i], 6); break }
        }
      }
      gsub(/"/, "", gid)
      sub(/^gene:/, "", gid)
      if (gid == "") {
        printf("ERROR: gene feature without ID/gene_id/Name at %s:%s-%s\n", chrom, $4, $5) > "/dev/stderr"
        exit 1
      }
      print chrom, start, end, gid, "0", strand
    }
  '
}

n=0
for fa in "${PROT_DIR}"/*.fa "${PROT_DIR}"/*.faa "${PROT_DIR}"/*.fasta; do
  [[ -e "${fa}" ]] || continue
  base=$(basename "${fa}")
  acc="${base%.*}"
  # strip .fa.gz-style double extension
  if [[ "${acc}" == *.fa || "${acc}" == *.faa || "${acc}" == *.fasta ]]; then
    acc="${acc%.*}"
  fi

  dest="${ANALYSIS}/accessions/${acc}"
  mkdir -p "${dest}"
  ln -sfn "${fa}" "${dest}/${acc}.fa"
  rm -f "${dest}/${acc}.bed"

  if [[ -n "${BED_DIR}" ]]; then
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
    to_genetribe_bed6 "${bed}" >"${dest}/${acc}.bed"
  else
    gff=""
    for cand in \
      "${GFF_DIR}/${acc}.gff3" "${GFF_DIR}/${acc}.gff" \
      "${GFF_DIR}/${acc}.gff3.gz" "${GFF_DIR}/${acc}.gff.gz"
    do
      if [[ -f "${cand}" ]]; then
        gff="${cand}"
        break
      fi
    done
    [[ -n "${gff}" ]] || {
      echo "ERROR: no GFF for accession '${acc}' under ${GFF_DIR}" >&2
      exit 1
    }
    gff_to_bed6 "${gff}" >"${dest}/${acc}.bed"
  fi

  ncols=$("${AWK[@]}" 'NF && $1 !~ /^#/ { print NF; exit }' "${dest}/${acc}.bed")
  [[ "${ncols}" == "6" ]] || {
    echo "ERROR: ${dest}/${acc}.bed has ${ncols:-0} columns after conversion (need 6)" >&2
    exit 1
  }
  n_genes=$("${AWK[@]}" 'NF && $1 !~ /^#/ { c++ } END { print c+0 }' "${dest}/${acc}.bed")
  [[ "${n_genes}" -gt 0 ]] || {
    echo "ERROR: ${dest}/${acc}.bed has 0 gene rows" >&2
    exit 1
  }

  cp -f "${CHRLIST}" "${dest}/${acc}.chrlist"
  n=$((n + 1))
  echo "  ${acc}: ${n_genes} genes → accessions/${acc}/"
done

[[ "${n}" -gt 0 ]] || {
  echo "ERROR: no protein FASTA files under ${PROT_DIR}" >&2
  exit 1
}
echo "Prepared ${n} accession(s) under ${ANALYSIS}/accessions/"
echo "Next: bash scripts/genetribe_check_ids.sh ${ANALYSIS}"
