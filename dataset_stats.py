"""
Summarize an OrthoFinder N0.tsv: accessions (columns), OGs, HOGs, gene rows.

Used by compute_dataset_stats.py and (after DB build) by build_index.py.
"""
from __future__ import annotations

import csv
import sys
from typing import Any


def summarize_n0_tsv(n0_path: str) -> dict[str, int]:
    """
    Count accessions from header, HOG rows, distinct orthogroups (OG), and genes
    (one per gene ID in cells, comma-separated lists allowed).
    """
    accessions = 0
    hogs = 0
    ogs: set[str] = set()
    genes = 0

    with open(n0_path, "r", newline="") as f:
        try:
            csv.field_size_limit(sys.maxsize)
        except OverflowError:
            csv.field_size_limit(2147483647)

        reader = csv.reader(f, delimiter="\t")
        header = next(reader)

        accession_cols: list[tuple[int, str]] = []
        for i, col in enumerate(header[3:], start=3):
            acc_name = col.replace(".primary.protein", "")
            accession_cols.append((i, acc_name))
        accessions = len(accession_cols)

        for row in reader:
            if len(row) < 3:
                continue
            og = row[1].strip()
            if og:
                ogs.add(og)
            hogs += 1
            for col_idx, _acc_name in accession_cols:
                if col_idx >= len(row):
                    continue
                cell = row[col_idx].strip()
                if not cell:
                    continue
                for gene_id in cell.split(", "):
                    if gene_id.strip():
                        genes += 1

    return {
        "accessions": accessions,
        "ogs": len(ogs),
        "hogs": hogs,
        "genes": genes,
    }


def write_dataset_stats_tsv(stats: dict[str, int], out_path: str) -> None:
    """Write a one-row TSV with header accessions, ogs, hogs, genes."""
    keys = ("accessions", "ogs", "hogs", "genes")
    with open(out_path, "w", newline="") as f:
        f.write("\t".join(keys) + "\n")
        f.write("\t".join(str(stats[k]) for k in keys) + "\n")


def stats_from_cursor(cur: Any) -> dict[str, int]:
    """Count rows using an open sqlite3 cursor (same connection as the build)."""
    cur.execute("SELECT COUNT(DISTINCT accession) FROM genes")
    accessions = int(cur.fetchone()[0])
    cur.execute("SELECT COUNT(DISTINCT og) FROM genes")
    ogs = int(cur.fetchone()[0])
    cur.execute("SELECT COUNT(DISTINCT hog) FROM hog_info")
    hogs = int(cur.fetchone()[0])
    cur.execute("SELECT COUNT(*) FROM genes")
    genes = int(cur.fetchone()[0])
    return {
        "accessions": accessions,
        "ogs": ogs,
        "hogs": hogs,
        "genes": genes,
    }


def stats_from_sqlite_db(db_path: str) -> dict[str, int] | None:
    """Read counts from a pan*.db built by build_index (after genes are loaded)."""
    import sqlite3

    try:
        conn = sqlite3.connect(db_path)
        cur = conn.cursor()
        out = stats_from_cursor(cur)
        conn.close()
    except (sqlite3.Error, OSError, TypeError, ValueError):
        return None
    return out
