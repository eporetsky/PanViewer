#!/usr/bin/env python3
"""
Parse an OrthoFinder N0 TSV and write dataset_stats.tsv (accessions, ogs, hogs, genes).

Example:
  python compute_dataset_stats.py --n0 input/wheat/N0.tsv --out input/wheat/dataset_stats.tsv
  python compute_dataset_stats.py --n0 input/barley/BPGv2_N0.tsv --out input/barley/dataset_stats.tsv
"""
import argparse
import os
import sys

from dataset_stats import summarize_n0_tsv, write_dataset_stats_tsv


def main() -> int:
    p = argparse.ArgumentParser(description="Summarize N0.tsv → dataset_stats.tsv")
    p.add_argument("--n0", required=True, help="Path to OrthoFinder N0.tsv")
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
        f"{stats['accessions']} accessions, {stats['ogs']} OGs, "
        f"{stats['hogs']:,} HOGs, {stats['genes']:,} genes"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
