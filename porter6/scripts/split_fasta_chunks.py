#!/usr/bin/env python3
"""Split a FASTA into multiple files with at most --chunk-size records each."""
from __future__ import annotations

import argparse
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, required=True, help="Input FASTA")
    ap.add_argument(
        "--chunk-size",
        type=int,
        default=100_000,
        help="Max sequences per output file (default: 100000)",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Directory for batch_000001.fasta, ...",
    )
    ap.add_argument(
        "--prefix",
        default="batch_",
        help="Output filename prefix (default: batch_)",
    )
    args = ap.parse_args()

    if args.chunk_size < 1:
        raise SystemExit("--chunk-size must be >= 1")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    batch_idx = 0
    n_in_batch = 0
    out_handle = None

    def open_next() -> None:
        nonlocal batch_idx, out_handle, n_in_batch
        if out_handle is not None:
            out_handle.close()
        batch_idx += 1
        path = args.out_dir / f"{args.prefix}{batch_idx:06d}.fasta"
        out_handle = open(path, "w", encoding="utf-8")
        n_in_batch = 0
        print(f"Writing {path.name}")

    open_next()

    for rec in SeqIO.parse(args.input, "fasta"):
        if n_in_batch >= args.chunk_size:
            open_next()
        rec.seq = Seq(str(rec.seq).replace("*", ""))
        SeqIO.write(rec, out_handle, "fasta")
        n_in_batch += 1

    if out_handle is not None:
        out_handle.close()

    print(f"Done: {batch_idx} file(s) in {args.out_dir}")


if __name__ == "__main__":
    main()
