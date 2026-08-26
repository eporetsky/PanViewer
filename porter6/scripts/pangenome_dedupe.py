#!/usr/bin/env python3
"""Deduplicate protein sequences across FASTA files; write unique FASTA + mapping TSV.

Mapping columns: orig_idx, source_file, header_id, uniq_id, seq_hash
  uniq_id is the FASTA header written for that sequence (Porter6 uses this as `id` in outputs).
"""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

from Bio import SeqIO


def normalize_seq(seq: str) -> str:
    s = str(seq).upper().replace(" ", "").replace("\n", "")
    return s


def seq_hash(seq: str) -> str:
    return hashlib.sha256(seq.encode("ascii", errors="replace")).hexdigest()[:16]


def iter_fasta_paths(root: Path) -> list[Path]:
    paths: list[Path] = []
    for pat in ("*.fa", "*.faa", "*.fasta", "*.FA", "*.FAA", "*.FASTA"):
        paths.extend(sorted(root.glob(pat)))
    return sorted(set(paths))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing .fa / .faa / .fasta files (non-recursive)",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Output directory",
    )
    ap.add_argument(
        "--unique-fasta",
        default="unique.fasta",
        help="Filename for deduplicated sequences (default: unique.fasta)",
    )
    ap.add_argument(
        "--map-tsv",
        default="orig_to_uniq.tsv",
        help="Filename for mapping TSV (default: orig_to_uniq.tsv)",
    )
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fasta_paths = iter_fasta_paths(args.input_dir)
    if not fasta_paths:
        raise SystemExit(f"No FASTA files found in {args.input_dir}")

    # seq_hash -> uniq_id (first time we see this sequence)
    hash_to_uniq: dict[str, str] = {}
    uniq_fasta_path = args.out_dir / args.unique_fasta
    map_path = args.out_dir / args.map_tsv

    n_in = 0
    n_unique = 0

    with open(map_path, "w", encoding="utf-8") as map_f:
        map_f.write("orig_idx\tsource_file\theader_id\tuniq_id\tseq_hash\n")
        with open(uniq_fasta_path, "w", encoding="utf-8") as ufa:
            for fpath in fasta_paths:
                for rec in SeqIO.parse(fpath, "fasta"):
                    seq = normalize_seq(rec.seq)
                    if not seq:
                        continue
                    n_in += 1
                    h = seq_hash(seq)
                    if h not in hash_to_uniq:
                        uniq_id = f"u{h}"
                        hash_to_uniq[h] = uniq_id
                        ufa.write(f">{uniq_id}\n")
                        for i in range(0, len(seq), 60):
                            ufa.write(seq[i : i + 60] + "\n")
                        n_unique += 1
                    else:
                        uniq_id = hash_to_uniq[h]
                    map_f.write(
                        f"{n_in}\t{fpath.name}\t{rec.id}\t{uniq_id}\t{h}\n"
                    )

    print(f"Input records: {n_in}")
    print(f"Unique sequences: {n_unique}")
    print(f"Wrote {uniq_fasta_path}")
    print(f"Wrote {map_path}")


if __name__ == "__main__":
    main()
