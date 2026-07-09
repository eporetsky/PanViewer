"""
Summarize index statistics: accessions, pan-gene clusters, gene rows.

Used by ``compute_dataset_stats.py`` (legacy tabular inputs) and ``build_index.py``.
"""
from __future__ import annotations

import csv
import os
import sys
from typing import Any

CORE_STATS_KEYS = ("accessions", "pan_genes", "genes")
OPTIONAL_TABLE_STATS_KEYS = (
    "gene_coords",
    "protein_seq_map",
    "cds_seq_map",
    "porter6_map",
)
STATS_KEYS = CORE_STATS_KEYS + OPTIONAL_TABLE_STATS_KEYS


def summarize_n0_tsv(n0_path: str) -> dict[str, int]:
    """
    Count accessions from header, cluster rows (one per TSV row after the header),
    and genes (one per gene ID in cells; comma-separated lists allowed).

    Legacy ``N0.tsv``-style pan-gene table layout; output uses the same keys as SQLite stats.
    """
    accessions = 0
    cluster_rows = 0
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
            cluster_rows += 1
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
        "pan_genes": cluster_rows,
        "genes": genes,
    }


def write_dataset_stats_tsv(stats: dict[str, int], out_path: str) -> None:
    """Write a one-row TSV with header accessions, pan_genes, genes."""
    keys = ("accessions", "pan_genes", "genes")
    with open(out_path, "w", newline="") as f:
        f.write("\t".join(keys) + "\n")
        f.write("\t".join(str(stats[k]) for k in keys) + "\n")


def write_combined_stats_tsv(
    out_path: str, rows: list[tuple[str, dict[str, int]]]
) -> None:
    """
    Multi-species stats for ``database/stats.tsv``:
    one row per species with core counts plus optional table row totals.
    """
    header = ("species",) + STATS_KEYS
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for species_id, stats in rows:
            sid = (species_id or "").strip().lower()
            line = [sid] + [str(int(stats[k])) for k in STATS_KEYS]
            f.write("\t".join(line) + "\n")


def _row_to_stats_dict(r: dict[str, str]) -> dict[str, int] | None:
    """Normalize a stats.tsv row to the full stats key set."""
    low = {((k or "").strip().lower()): (v or "").strip() for k, v in r.items()}
    acc_s = low.get("accessions") or ""
    genes_s = low.get("genes") or ""
    pan_s = (low.get("pan_genes") or "").strip()
    if not all((acc_s, pan_s, genes_s)):
        return None
    try:
        out = {
            "accessions": int(acc_s),
            "pan_genes": int(pan_s),
            "genes": int(genes_s),
        }
        for key in OPTIONAL_TABLE_STATS_KEYS:
            raw = low.get(key) or ""
            out[key] = int(raw) if raw else 0
        return out
    except ValueError:
        return None


def merge_write_combined_stats_tsv(
    out_path: str, updates: list[tuple[str, dict[str, int]]]
) -> None:
    """Merge into ``database/stats.tsv``, keeping species rows not rebuilt this run."""
    merged: dict[str, dict[str, int]] = {}
    if os.path.isfile(out_path):
        try:
            with open(out_path, newline="", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                fnames = [((x or "").strip().lower()) for x in (reader.fieldnames or [])]
                if "species" in fnames:
                    for row in reader:
                        r = {((k or "").strip().lower()): v for k, v in row.items()}
                        sid = (r.get("species") or "").strip().lower()
                        if not sid:
                            continue
                        parsed = _row_to_stats_dict(r)
                        if parsed:
                            merged[sid] = parsed
                elif "accessions" in fnames and "genes" in fnames:
                    row = next(reader, None)
                    if row:
                        r = {((k or "").strip().lower()): v for k, v in row.items()}
                        parsed = _row_to_stats_dict(r)
                        if parsed:
                            merged["wheat"] = parsed
        except (OSError, ValueError, TypeError):
            merged = {}
    for species_id, stats in updates:
        sid = (species_id or "").strip().lower()
        if sid:
            merged[sid] = {k: int(stats[k]) for k in STATS_KEYS}
    out_rows = [(s, merged[s]) for s in sorted(merged.keys())]
    write_combined_stats_tsv(out_path, out_rows)


def _table_exists(cur: Any, table: str) -> bool:
    cur.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=? LIMIT 1",
        (table,),
    )
    return cur.fetchone() is not None


def _table_row_count(cur: Any, table: str) -> int:
    if not _table_exists(cur, table):
        return 0
    cur.execute(f"SELECT COUNT(*) FROM {table}")
    return int(cur.fetchone()[0])


def _genes_pangene_column(cur: Any) -> str | None:
    cur.execute("PRAGMA table_info(genes)")
    cols = {str(row[1]) for row in cur.fetchall()}
    if "pangene" in cols:
        return "pangene"
    if "pan_gene" in cols:
        return "pan_gene"
    return None


def stats_from_cursor(cur: Any) -> dict[str, int]:
    """Count rows using an open sqlite3 cursor (same connection as the build)."""
    cur.execute("SELECT COUNT(DISTINCT accession) FROM genes")
    accessions = int(cur.fetchone()[0])
    if cur.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='pangene_info' LIMIT 1"
    ).fetchone():
        cur.execute("SELECT COUNT(*) FROM pangene_info")
    else:
        pan_col = _genes_pangene_column(cur)
        if not pan_col:
            raise ValueError("genes table has no pangene or pan_gene column")
        cur.execute(f"SELECT COUNT(DISTINCT {pan_col}) FROM genes")
    pan_genes = int(cur.fetchone()[0])
    cur.execute("SELECT COUNT(*) FROM genes")
    genes = int(cur.fetchone()[0])
    out = {
        "accessions": accessions,
        "pan_genes": pan_genes,
        "genes": genes,
    }
    for table in OPTIONAL_TABLE_STATS_KEYS:
        out[table] = _table_row_count(cur, table)
    return out


def stats_from_sqlite_db(db_path: str) -> dict[str, int] | None:
    """Read counts from a species ``*.db`` built by ``build_index.py``."""
    import sqlite3

    try:
        conn = sqlite3.connect(db_path)
        cur = conn.cursor()
        out = stats_from_cursor(cur)
        conn.close()
    except (sqlite3.Error, OSError, TypeError, ValueError):
        return None
    return out


def discover_species_db_paths(db_dir: str) -> list[tuple[str, str]]:
    """Return ``(species_id, db_path)`` for each ``<species>.db`` under ``db_dir``."""
    if not os.path.isdir(db_dir):
        return []
    out: list[tuple[str, str]] = []
    for fn in sorted(os.listdir(db_dir)):
        if not fn.endswith(".db"):
            continue
        stem = os.path.splitext(fn)[0].strip().lower()
        if not stem or stem.startswith("."):
            continue
        out.append((stem, os.path.join(db_dir, fn)))
    return out


def collect_stats_from_database_dir(db_dir: str) -> list[tuple[str, dict[str, int]]]:
    """Scan every species ``*.db`` and return stats rows (skips unreadable DBs)."""
    rows: list[tuple[str, dict[str, int]]] = []
    for species_id, db_path in discover_species_db_paths(db_dir):
        stats = stats_from_sqlite_db(db_path)
        if stats is not None:
            rows.append((species_id, stats))
    return rows
