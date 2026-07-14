#!/usr/bin/env bash
# Append .$fa to ingest outputs under 02_fasta_nuc so globs like 02_fasta_nuc/*.$fa work.
# Safe to run repeatedly (only touches lines that still look unpatched).
#
# Run from the repo root (same cwd as the Makefile).

set -euo pipefail

PAN_SH="$(pwd)/pandagma-upstream/bin/pandagma-pan.sh"

if [[ ! -f "$PAN_SH" ]]; then
  echo "ERROR: not found: $PAN_SH — cd to repo root first." >&2
  exit 1
fi

echo "Patching: $PAN_SH"

perl -i.bak -pe '
  s/^(\s*-out 02_fasta_nuc\/"\$file_base")$/\1."\$fa"/;
  s/(basename 02_fasta_nuc\/"\$file_base")\s*\|/\1."\$fa" \|/;
  s/(mmseqs easy-search "02_fasta_nuc\/\$file_base") /\1."\$fa" /;
  s|echo "Extra: 02_fasta_nuc/\$file_base"|echo "Extra: 02_fasta_nuc/\$file_base.\$fa"|;
' "$PAN_SH"

rm -f "${PAN_SH}.bak"
echo "Done."
