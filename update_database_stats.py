#!/usr/bin/env python3
"""
Refresh ``database/stats.tsv`` from existing ``database/<species>.db`` files.

Use after copying or rebuilding databases without re-running the full
``build_index.py`` pipeline, or to backfill extended table counts.

  python update_database_stats.py
  python update_database_stats.py --database-dir /path/to/database
"""
from __future__ import annotations

import argparse
import os
import sys

from dataset_stats import (
    discover_species_db_paths,
    stats_from_sqlite_db,
    write_combined_stats_tsv,
)


def _format_stats_line(species_id: str, stats: dict[str, int]) -> str:
    return (
        f"  {species_id}: {stats['accessions']} accessions, "
        f"{stats['pan_genes']:,} pan-genes, {stats['genes']:,} genes, "
        f"{stats['gene_coords']:,} gene_coords, "
        f"{stats['protein_seq_map']:,} protein_seq_map, "
        f"{stats['cds_seq_map']:,} cds_seq_map, "
        f"{stats['porter6_map']:,} porter6_map"
    )


def main(argv: list[str] | None = None) -> int:
    repo_root = os.path.dirname(os.path.abspath(__file__))
    default_db_dir = os.path.join(repo_root, "database")
    default_out = os.path.join(default_db_dir, "stats.tsv")

    p = argparse.ArgumentParser(
        description="Parse database/*.db and write an updated database/stats.tsv"
    )
    p.add_argument(
        "--database-dir",
        default=default_db_dir,
        help=f"Directory containing <species>.db files (default: {default_db_dir})",
    )
    p.add_argument(
        "--out",
        default=default_out,
        help=f"Output TSV path (default: {default_out})",
    )
    args = p.parse_args(argv)

    db_dir = os.path.abspath(args.database_dir)
    out_path = os.path.abspath(args.out)

    dbs = discover_species_db_paths(db_dir)
    if not dbs:
        print(f"No *.db files under {db_dir}", file=sys.stderr)
        return 1

    rows: list[tuple[str, dict[str, int]]] = []
    failed: list[str] = []
    for species_id, db_path in dbs:
        stats = stats_from_sqlite_db(db_path)
        if stats is None:
            failed.append(species_id)
        else:
            rows.append((species_id, stats))

    if not rows:
        print("No readable databases; stats.tsv not updated.", file=sys.stderr)
        return 1

    if failed:
        print(
            f"Warning: could not read stats from: {', '.join(failed)}",
            file=sys.stderr,
        )

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    write_combined_stats_tsv(out_path, rows)

    print(f"Updated {out_path} ({len(rows)} species)")
    for species_id, stats in rows:
        print(_format_stats_line(species_id, stats))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
