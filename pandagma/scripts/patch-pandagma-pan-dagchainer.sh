#!/usr/bin/env bash
# Patch pandagma-pan.sh run_dagchainer for reliable Slurm/batch runs:
#   - Limit parallel DAG jobs (DAG_NPROC, default 8) — 48× dagchainer often OOMs → empty .aligncoords
#   - Log each pair instead of discarding stdout (1>/dev/null)
#   - Skip pairs that already have non-empty .aligncoords (resume-friendly)
#   - Warn when output is still empty after a run
#
# Run from repo root after make setup (targets pandagma-upstream/bin/pandagma-pan.sh).

set -euo pipefail

PAN_SH="$(pwd)/pandagma-upstream/bin/pandagma-pan.sh"

if [[ ! -f "$PAN_SH" ]]; then
  echo "ERROR: not found: $PAN_SH — cd to repo root first." >&2
  exit 1
fi

if grep -q 'DAG_NPROC' "$PAN_SH" 2>/dev/null; then
  echo "Already patched (DAG_NPROC present): $PAN_SH"
  exit 0
fi

cp -a "$PAN_SH" "${PAN_SH}.pre-dagchainer-patch"

python3 <<'PY'
from pathlib import Path
import re

path = Path("pandagma-upstream/bin/pandagma-pan.sh")
text = path.read_text()
old = r'''run_dagchainer\(\) \{
  # Identify syntenic blocks, using DAGchainer
  cd "\$\{WORK_DIR\}" \|\| exit
  dagchainer_args='-M 50 -E 1e-5 -A 6 -s'  # -g and -D are calculated from the data
  echo; echo "Run DAGchainer using arguments \"\$\{dagchainer_args\}\" \(-g and -D are calculated from the data\)"
  # Check and preemptively remove malformed \\\*_matches\.file, which can result from an aborted run
  if \[ -f 04_dag/\\\*_matches\.tsv \]; then rm 04_dag/\\\*_matches\.tsv; fi
  for match_path in 04_dag/\*_matches\.tsv; do
    align_file=\$\(basename "\$match_path" _matches\.tsv\)
    qryfile=\$\(echo "\$align_file" \| perl -pe 's/\(\S\+\)\.x\..+/\$1/'\)
    sbjfile=\$\(echo "\$align_file" \| perl -pe 's/\S\+\.x\.\(\S\+\)/\$1/'\)

    echo "Find average distance between genes for the query and subject files: "
    echo "  \$qryfile and \$sbjfile"
    ave_gene_gap=\$\(cat 02_fasta_nuc/"\$qryfile"\."\$fa" 02_fasta_nuc/"\$sbjfile"\."\$fa" \| 
                     awk '\$1~/\^>/ \{print substr\(\$1,2\)\}' \| perl -pe 's/__/\t/g' \| sort -k1,1 -k3n,3n \|
                     awk '\$1 == prev1 && \$3 > prev4 \{sum\+\=\$3-prev4; ct\+\+; prev1=\$1; prev3=\$3; prev4=\$4\};
                          \$1 != prev1 \|\| \$3 <= prev4 \{prev1=\$1; prev3=\$3; prev4=\$4\}; 
                          END\{print 100\*int\(sum/ct/100\)\}'\)
    max_gene_gap=\$\(\( ave_gene_gap \* 20 \)\)

    echo "Running DAGchainer on comparison: \$align_file"
    echo "  Calculated DAGchainer parameters: -g \(ave_gene_gap\): \$ave_gene_gap -D \(max_gene_gap\): \$max_gene_gap"; echo
    echo " run_DAG_chainer\.pl \$dagchainer_args  -g \$ave_gene_gap -D \$max_gene_gap -i \"\$\{OLDPWD\}/\$\{match_path\}\""

    # run_DAG_chainer\.pl writes temp files to cwd;
    # use per-process temp directory to avoid any data race
    \(
      tmpdir=\$\(mktemp -d\)
      cd "\$\{tmpdir\}" \|\| exit
      run_DAG_chainer\.pl "\$dagchainer_args"  -g "\$ave_gene_gap" -D "\$max_gene_gap" -i "\$\{OLDPWD\}/\$\{match_path\}" 1>/dev/null
      rmdir "\$\{tmpdir\}"
    \) &
    # allow to execute up to \$NPROC in parallel
    \[ "\$\(jobs -r -p \| wc -l\)" -ge "\$\{NPROC\}" \] && wait -n
  done
  wait # wait for last jobs to finish

  if \[ "\$strict_synt" -eq 1 \]; then'''

# Simpler: find function start/end by line numbers from known file
lines = text.splitlines(keepends=True)
start = end = None
for i, line in enumerate(lines):
    if line.startswith("run_dagchainer() {"):
        start = i
    if start is not None and line.startswith("run_mcl() {"):
        end = i
        break
if start is None or end is None:
    raise SystemExit("Could not locate run_dagchainer() in pandagma-pan.sh")

new_func = r'''run_dagchainer() {
  # Identify syntenic blocks, using DAGchainer
  cd "${WORK_DIR}" || exit
  mkdir -p stats/dagchainer_logs
  : "${DAG_NPROC:=8}"
  : "${DAG_SKIP_EXISTING:=1}"
  echo; echo "Run DAGchainer (parallel cap DAG_NPROC=${DAG_NPROC}, skip existing=${DAG_SKIP_EXISTING})"
  dag_A="${dag_A:-6}"
  dag_M="${dag_M:-50}"
  dag_E="${dag_E:-1e-5}"
  dag_gap_mult="${dag_gap_mult:-20}"
  dagchainer_args="-M ${dag_M} -E ${dag_E} -A ${dag_A} -s"
  echo "  dagchainer_args: ${dagchainer_args} (-g and -D per pair; gap mult ${dag_gap_mult})"
  if [ -f 04_dag/\*_matches.tsv ]; then rm 04_dag/\*_matches.tsv; fi
  for match_path in 04_dag/*_matches.tsv; do
    align_file=$(basename "$match_path" _matches.tsv)
    align_out="${match_path}.aligncoords"
    qryfile=$(echo "$align_file" | perl -pe 's/(\S+)\.x\..+/$1/')
    sbjfile=$(echo "$align_file" | perl -pe 's/\S+\.x\.(\S+)/$1/')

    if [[ "${DAG_SKIP_EXISTING}" == 1 && -f "${align_out}" ]]; then
      ac_lines=$(grep -cve '^#' -e '^$' "${align_out}" 2>/dev/null || echo 0)
      if [[ "${ac_lines}" -gt 0 ]]; then
        echo "  Skip ${align_file}: ${align_out} already has ${ac_lines} data lines"
        continue
      fi
    fi

    echo "Find average distance between genes for the query and subject files: "
    echo "  $qryfile and $sbjfile"
    ave_gene_gap=$(cat 02_fasta_nuc/"$qryfile"."$fa" 02_fasta_nuc/"$sbjfile"."$fa" |
                     awk '$1~/^>/ {print substr($1,2)}' | perl -pe 's/__/\t/g' | sort -k1,1 -k3n,3n |
                     awk '$1 == prev1 && $3 > prev4 {sum+=$3-prev4; ct++; prev1=$1; prev3=$3; prev4=$4};
                          $1 != prev1 || $3 <= prev4 {prev1=$1; prev3=$3; prev4=$4};
                          END{if (ct>0) print 100*int(sum/ct/100); else print 0}')
    max_gene_gap=$(( ave_gene_gap * dag_gap_mult ))

    echo "Running DAGchainer on comparison: $align_file"
    echo "  Calculated DAGchainer parameters: -g (ave_gene_gap): $ave_gene_gap -D (max_gene_gap): $max_gene_gap"

    (
      tmpdir=$(mktemp -d)
      cd "${tmpdir}" || exit
      log="${OLDPWD}/stats/dagchainer_logs/${align_file}.log"
      if ! run_DAG_chainer.pl ${dagchainer_args} -g "$ave_gene_gap" -D "$max_gene_gap" \
          -i "${OLDPWD}/${match_path}" >>"${log}" 2>&1; then
        echo "WARNING: run_DAG_chainer.pl exited non-zero for ${align_file} (see ${log})" >&2
      fi
      if [[ ! -s "${OLDPWD}/${align_out}" ]]; then
        shopt -s nullglob
        for stray in ./*.aligncoords ./*aligncoords; do
          if [[ -f "${stray}" ]]; then
            mv "${stray}" "${OLDPWD}/${align_out}"
            break
          fi
        done
        shopt -u nullglob
      fi
      rm -rf "${tmpdir}"
      ac_lines=$(grep -cve '^#' -e '^$' "${OLDPWD}/${align_out}" 2>/dev/null || echo 0)
      if [[ "${ac_lines}" -eq 0 ]]; then
        echo "WARNING: empty aligncoords for ${align_file} (log: ${log})" >&2
      fi
    ) &
    if [[ $(jobs -r -p | wc -l) -ge ${DAG_NPROC} ]]; then wait -n; fi
  done
  wait

  if [ "$strict_synt" -eq 1 ]; then
    cat 04_dag/*.aligncoords | awk '$1!~/^#/ {print $2 "\t" $6}' |
      awk 'NF==2' | sort -u > 05_filtered_pairs.tsv
  else
    cat 04_dag/*_matches.tsv 04_dag/*.aligncoords | awk '$1!~/^#/ {print $2 "\t" $6}' |
      awk 'NF==2' | sort -u > 05_filtered_pairs.tsv
  fi
}

##########
'''

lines[start:end] = [new_func]
path.write_text("".join(lines))
print(f"Patched run_dagchainer() in {path} (lines {start+1}-{end})")
PY

echo "Done. Backup: ${PAN_SH}.pre-dagchainer-patch"
echo "Re-run: make patch-pandagma-dagchainer  (after git pull of upstream)"
