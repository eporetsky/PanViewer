# Shared DAGchainer helpers (source from shard scripts; do not execute directly).
# Expects: WORK_DIR, CONF, DATA_DIR set; caller cd's to WORK_DIR before run_dag_pairs.

pandagma_dagchainer_load_config() {
  if [[ ! -d "${WORK_DIR}/02_fasta_nuc" ]]; then
    echo "ERROR: ${WORK_DIR}/02_fasta_nuc missing (run ingest first)." >&2
    return 1
  fi
  if ! command -v run_DAG_chainer.pl >/dev/null 2>&1; then
    echo "ERROR: run_DAG_chainer.pl not on PATH (conda activate pandagma?)." >&2
    return 1
  fi

  : "${DAG_NPROC:=4}"
  : "${DAG_SKIP_EXISTING:=1}"

  declare -g -a _dag_cds_files
  cd "${DATA_DIR}" || return 1
  # shellcheck disable=SC1090
  source "${CONF}"

  dag_A="${dag_A:-6}"
  dag_M="${dag_M:-50}"
  dag_E="${dag_E:-1e-5}"
  dag_gap_mult="${dag_gap_mult:-20}"
  dagchainer_args="-M ${dag_M} -E ${dag_E} -A ${dag_A} -s"

  mapfile -t _dag_cds_files < <(realpath --canonicalize-existing "${cds_files[@]}")
  local fasta_file
  fasta_file=$(basename "${_dag_cds_files[0]}" .gz)
  fa="${fasta_file##*.}"

  mkdir -p "${WORK_DIR}/stats/dagchainer_logs" "${WORK_DIR}/04_dag"
  cd "${WORK_DIR}" || return 1
}

pandagma_dagchainer_run_one_pair() {
  local match_path=$1
  local align_file align_out qryfile sbjfile ave_gene_gap max_gene_gap log ac_lines

  align_file=$(basename "$match_path" _matches.tsv)
  align_out="${match_path}.aligncoords"

  if [[ "${DAG_SKIP_EXISTING}" == 1 && -f "${align_out}" ]]; then
    ac_lines=$(grep -cve '^#' -e '^$' "${align_out}" 2>/dev/null || echo 0)
    if [[ "${ac_lines}" -gt 0 ]]; then
      echo "  Skip ${align_file} (${ac_lines} lines)"
      return 0
    fi
  fi

  qryfile=$(echo "$align_file" | perl -pe 's/(\S+)\.x\..+/$1/')
  sbjfile=$(echo "$align_file" | perl -pe 's/\S+\.x\.(\S+)/$1/')

  ave_gene_gap=$(cat "02_fasta_nuc/${qryfile}.${fa}" "02_fasta_nuc/${sbjfile}.${fa}" |
    awk '$1~/^>/ {print substr($1,2)}' | perl -pe 's/__/\t/g' | sort -k1,1 -k3n,3n |
    awk '$1 == prev1 && $3 > prev4 {sum+=$3-prev4; ct++; prev1=$1; prev3=$3; prev4=$4};
         $1 != prev1 || $3 <= prev4 {prev1=$1; prev3=$3; prev4=$4};
         END{if (ct>0) print 100*int(sum/ct/100); else print 0}')
  max_gene_gap=$(( ave_gene_gap * dag_gap_mult ))

  echo "  Run ${align_file} -g ${ave_gene_gap} -D ${max_gene_gap}"
  log="stats/dagchainer_logs/${align_file}.log"
  local tmpdir
  tmpdir=$(mktemp -d)
  (
    cd "${tmpdir}" || exit 1
    if ! run_DAG_chainer.pl ${dagchainer_args} -g "${ave_gene_gap}" -D "${max_gene_gap}" \
        -i "${OLDPWD}/${match_path}" >>"${OLDPWD}/${log}" 2>&1; then
      echo "WARNING: run_DAG_chainer.pl failed for ${align_file} (see ${log})" >&2
    fi
    if [[ ! -s "${OLDPWD}/${align_out}" ]]; then
      shopt -s nullglob
      for stray in ./*.aligncoords ./*aligncoords; do
        [[ -f "${stray}" ]] || continue
        mv "${stray}" "${OLDPWD}/${align_out}"
        break
      done
      shopt -u nullglob
    fi
    rm -rf "${tmpdir}"
  )
  ac_lines=$(grep -cve '^#' -e '^$' "${align_out}" 2>/dev/null || echo 0)
  if [[ "${ac_lines}" -eq 0 ]]; then
    echo "WARNING: empty aligncoords for ${align_file} (see ${log})" >&2
  else
    echo "  Done ${align_file}: ${ac_lines} data lines"
  fi
}

pandagma_dagchainer_run_pairs() {
  local -a matches=("$@")
  local match_path ok=0 empty=0

  for match_path in "${matches[@]}"; do
    pandagma_dagchainer_run_one_pair "${match_path}" &
    while [[ $(jobs -r -p | wc -l) -ge ${DAG_NPROC} ]]; do
      wait -n
    done
  done
  wait

  for match_path in "${matches[@]}"; do
    local ac="${match_path}.aligncoords"
    if [[ -f "${ac}" ]] && [[ $(grep -cve '^#' -e '^$' "${ac}" 2>/dev/null || echo 0) -gt 0 ]]; then
      ok=$((ok + 1))
    else
      empty=$((empty + 1))
    fi
  done
  echo "Batch done: ${ok} with synteny, ${empty} empty/missing (of ${#matches[@]} pairs)."
}
