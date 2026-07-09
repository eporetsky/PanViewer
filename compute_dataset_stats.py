#!/usr/bin/env python3
"""
Parse a legacy ``N0.tsv`` (pan-gene / gene matrix layout) and write ``dataset_stats.tsv``
(accessions, pan_genes, genes).

Example:
  python compute_dataset_stats.py --n0 input/wheat/N0.tsv --out input/wheat/dataset_stats.tsv
"""
import argparse
import os
import sys

from dataset_stats import summarize_n0_tsv, write_dataset_stats_tsv


def main() -> int:
    p = argparse.ArgumentParser(description="Summarize N0.tsv → dataset_stats.tsv (legacy layout)")
    p.add_argument("--n0", required=True, help="Path to N0.tsv (legacy pan-gene matrix layout)")
    p.add_argument(
        "--out",
        required=True,
        help="Output TSV path (e.g. input/wheat/dataset_stats.tsv)",
    )
    args = p.parse_args()
    n0 = os.path.abspath(args.n0)
    out = os.path.abspath(args.out)
    if not os.path.isfile(n0):
        print(f"Not found: {n0}", file=sys.stderr)
        return 1
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    stats = summarize_n0_tsv(n0)
    write_dataset_stats_tsv(stats, out)
    print(
        f"Wrote {out}: "
        f"{stats['accessions']} accessions, {stats['pan_genes']:,} pan-genes, "
        f"{stats['genes']:,} genes"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
