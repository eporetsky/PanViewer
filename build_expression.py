#!/usr/bin/env python3
"""
Build ``expression/expression.db`` from per-dataset folders under ``expression/``.

Each dataset directory (e.g. ``expression/wheat/``) may contain:

  meta.tsv / meta.tsv.gz     required — one row per sample
  cpm.tsv / cpm.tsv.gz       required — gene × sample CPM matrix
  deg.tsv / deg.tsv.gz       optional — long DEG table
  dataset.json               optional — ``genome``, ``label`` overrides

All rows present in these files are ingested (no sample/gene filtering). Prep
your TSVs upstream if you want a smaller database.

Runtime: if ``expression/expression.db`` exists, PanViewer uses it instead of the
PlantApp API (see ``expression_local.py`` / ``plantapp_omics.py``).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import sqlite3
import sys
import time
from collections import Counter
from typing import Any, Iterable, TextIO

from expression_local import (
    EXPRESSION_DB_NAME,
    SCHEMA_VERSION,
    encode_cpm_values,
    expression_db_path,
)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
EXPRESSION_DIR = os.path.join(BASE_DIR, "expression")

META_REQUIRED = ("sample_acc",)
META_OPTIONAL = (
    "experiment_acc",
    "organ",
    "group",
    "genome",
    "species",
    "stage",
    "inbred",
    "short_title",
    "title",
)
# Accept PlantApp typo / alias headers when reading sources.
META_ALIASES = {
    "sample_accession": "sample_acc",
    "experiment_accession": "experiment_acc",
    "Sample": "sample_acc",
    "sample": "sample_acc",
}

DEG_REQUIRED = ("gene_id",)
DEG_OPTIONAL = (
    "experiment_acc",
    "category",
    "comparison",
    "group1",
    "group2",
    "logFC",
    "FDR",
    "short_title",
    "title",
)
DEG_ALIASES = {
    "GeneID": "gene_id",
    "geneID": "gene_id",
    "GeneId": "gene_id",
    "logfc": "logFC",
    "LogFC": "logFC",
    "fdr": "FDR",
    "adj.P.Val": "FDR",
}


def _open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def _find_file(dataset_dir: str, stems: Iterable[str]) -> str | None:
    for stem in stems:
        for name in (f"{stem}.tsv", f"{stem}.tsv.gz", f"{stem}.csv", f"{stem}.csv.gz"):
            path = os.path.join(dataset_dir, name)
            if os.path.isfile(path):
                return path
    return None


def _normalize_header(name: str, aliases: dict[str, str]) -> str:
    n = (name or "").strip()
    if n in aliases:
        return aliases[n]
    return n


def _dialect_for(path: str) -> str:
    return "|" if path.endswith((".csv", ".csv.gz")) else "\t"


def _load_dataset_json(dataset_dir: str) -> dict[str, Any]:
    path = os.path.join(dataset_dir, "dataset.json")
    if not os.path.isfile(path):
        return {}
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise RuntimeError(f"{path}: expected a JSON object")
    return data


def _read_meta(path: str) -> list[dict[str, str]]:
    delim = "," if path.endswith((".csv", ".csv.gz")) else "\t"
    with _open_text(path) as f:
        reader = csv.DictReader(f, delimiter=delim)
        if not reader.fieldnames:
            raise RuntimeError(f"{path}: empty header")
        field_map = {
            raw: _normalize_header(raw, META_ALIASES) for raw in reader.fieldnames
        }
        rows: list[dict[str, str]] = []
        for i, raw in enumerate(reader, start=2):
            row = {field_map[k]: (v or "").strip() for k, v in raw.items() if k is not None}
            if not row.get("sample_acc"):
                raise RuntimeError(f"{path}:{i}: missing sample_acc")
            rows.append(row)
    if not rows:
        raise RuntimeError(f"{path}: no sample rows")
    return rows


def _majority_genome(rows: list[dict[str, str]]) -> str:
    counts = Counter((r.get("genome") or "").strip() for r in rows if (r.get("genome") or "").strip())
    if not counts:
        return ""
    return counts.most_common(1)[0][0]


def _ensure_schema(con: sqlite3.Connection) -> None:
    con.executescript(
        """
        CREATE TABLE meta (
            key TEXT PRIMARY KEY,
            value TEXT NOT NULL
        );
        CREATE TABLE datasets (
            dataset_id TEXT PRIMARY KEY,
            genome TEXT NOT NULL,
            label TEXT NOT NULL
        );
        CREATE TABLE samples (
            dataset_id TEXT NOT NULL,
            sample_idx INTEGER NOT NULL,
            sample_acc TEXT NOT NULL,
            experiment_acc TEXT,
            organ TEXT,
            "group" TEXT,
            genome TEXT,
            species TEXT,
            stage TEXT,
            inbred TEXT,
            short_title TEXT,
            title TEXT,
            PRIMARY KEY (dataset_id, sample_idx),
            UNIQUE (dataset_id, sample_acc)
        );
        CREATE TABLE experiments (
            dataset_id TEXT NOT NULL,
            experiment_acc TEXT NOT NULL,
            short_title TEXT,
            title TEXT,
            PRIMARY KEY (dataset_id, experiment_acc)
        );
        CREATE TABLE cpm (
            dataset_id TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            cpm_blob BLOB NOT NULL,
            PRIMARY KEY (dataset_id, gene_id)
        );
        CREATE TABLE deg (
            dataset_id TEXT NOT NULL,
            gene_id TEXT NOT NULL,
            experiment_acc TEXT NOT NULL,
            category TEXT,
            comparison TEXT,
            group1 TEXT,
            group2 TEXT,
            logFC REAL,
            FDR REAL
        );
        CREATE INDEX idx_cpm_gene ON cpm(gene_id);
        CREATE INDEX idx_deg_gene ON deg(dataset_id, gene_id);
        """
    )


def _ingest_cpm(
    con: sqlite3.Connection,
    dataset_id: str,
    path: str,
    sample_accs: list[str],
) -> int:
    delim = "," if path.endswith((".csv", ".csv.gz")) else "\t"
    sample_set = set(sample_accs)
    n_genes = 0
    batch: list[tuple[str, str, bytes]] = []

    with _open_text(path) as f:
        header = f.readline()
        if not header:
            raise RuntimeError(f"{path}: empty file")
        cols = header.rstrip("\n\r").split(delim)
        if len(cols) < 2:
            raise RuntimeError(f"{path}: expected gene_id + sample columns")
        gene_col = cols[0].strip()
        if gene_col.lower() not in ("gene_id", "geneid", "gene"):
            # Still accept, but warn via stderr
            print(
                f"  note: CPM first column is {gene_col!r} (expected gene_id)",
                file=sys.stderr,
            )
        col_sample = [c.strip() for c in cols[1:]]
        missing = [s for s in sample_accs if s not in col_sample]
        extra = [c for c in col_sample if c not in sample_set]
        if missing:
            raise RuntimeError(
                f"{path}: CPM matrix missing {len(missing)} sample(s) from meta "
                f"(e.g. {missing[:3]})"
            )
        if extra:
            print(
                f"  note: CPM has {len(extra)} sample column(s) not in meta; ignored",
                file=sys.stderr,
            )
        index_of = {acc: i for i, acc in enumerate(col_sample)}
        order_idx = [index_of[s] for s in sample_accs]

        for line_no, line in enumerate(f, start=2):
            line = line.rstrip("\n\r")
            if not line:
                continue
            parts = line.split(delim)
            if len(parts) != len(cols):
                raise RuntimeError(
                    f"{path}:{line_no}: expected {len(cols)} fields, got {len(parts)}"
                )
            gene_id = parts[0].strip()
            if not gene_id:
                raise RuntimeError(f"{path}:{line_no}: empty gene_id")
            values: list[float | None] = []
            for j in order_idx:
                cell = parts[j + 1].strip()
                if cell == "" or cell.upper() in ("NA", "NAN", "NULL", "."):
                    values.append(None)
                else:
                    try:
                        values.append(float(cell))
                    except ValueError as e:
                        raise RuntimeError(
                            f"{path}:{line_no}: bad CPM value {cell!r}"
                        ) from e
            batch.append((dataset_id, gene_id, encode_cpm_values(values)))
            n_genes += 1
            if len(batch) >= 500:
                con.executemany(
                    "INSERT INTO cpm(dataset_id, gene_id, cpm_blob) VALUES (?, ?, ?)",
                    batch,
                )
                batch.clear()

    if batch:
        con.executemany(
            "INSERT INTO cpm(dataset_id, gene_id, cpm_blob) VALUES (?, ?, ?)",
            batch,
        )
    return n_genes


def _parse_float(cell: str, path: str, line_no: int, field: str) -> float | None:
    cell = (cell or "").strip()
    if cell == "" or cell.upper() in ("NA", "NAN", "NULL", "."):
        return None
    try:
        return float(cell)
    except ValueError as e:
        raise RuntimeError(f"{path}:{line_no}: bad {field} value {cell!r}") from e


def _ingest_deg(con: sqlite3.Connection, dataset_id: str, path: str) -> tuple[int, int]:
    delim = "," if path.endswith((".csv", ".csv.gz")) else "\t"
    n_rows = 0
    experiments: dict[str, tuple[str, str]] = {}
    batch: list[tuple[Any, ...]] = []

    with _open_text(path) as f:
        reader = csv.DictReader(f, delimiter=delim)
        if not reader.fieldnames:
            raise RuntimeError(f"{path}: empty DEG header")
        field_map = {
            raw: _normalize_header(raw, DEG_ALIASES) for raw in reader.fieldnames
        }
        norm_fields = set(field_map.values())
        if "gene_id" not in norm_fields:
            raise RuntimeError(f"{path}: DEG table needs a gene_id / GeneID column")

        for line_no, raw in enumerate(reader, start=2):
            row = {field_map[k]: (v or "").strip() for k, v in raw.items() if k is not None}
            gene_id = row.get("gene_id") or ""
            if not gene_id:
                raise RuntimeError(f"{path}:{line_no}: empty gene_id")
            ea = row.get("experiment_acc") or ""
            if not ea:
                raise RuntimeError(f"{path}:{line_no}: empty experiment_acc")
            st = row.get("short_title") or ""
            title = row.get("title") or ""
            prev = experiments.get(ea)
            if prev is None:
                experiments[ea] = (st, title)
            elif (st or title) and (not prev[0] and not prev[1]):
                experiments[ea] = (st, title)

            log_fc = _parse_float(row.get("logFC") or "", path, line_no, "logFC")
            fdr = _parse_float(row.get("FDR") or "", path, line_no, "FDR")
            batch.append(
                (
                    dataset_id,
                    gene_id,
                    ea,
                    row.get("category") or "",
                    row.get("comparison") or "",
                    row.get("group1") or "",
                    row.get("group2") or "",
                    log_fc,
                    fdr,
                )
            )
            n_rows += 1
            if len(batch) >= 2000:
                con.executemany(
                    """
                    INSERT INTO deg(
                        dataset_id, gene_id, experiment_acc, category, comparison,
                        group1, group2, logFC, FDR
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    batch,
                )
                batch.clear()

    if batch:
        con.executemany(
            """
            INSERT INTO deg(
                dataset_id, gene_id, experiment_acc, category, comparison,
                group1, group2, logFC, FDR
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            batch,
        )

    con.executemany(
        """
        INSERT INTO experiments(dataset_id, experiment_acc, short_title, title)
        VALUES (?, ?, ?, ?)
        """,
        [(dataset_id, ea, st, title) for ea, (st, title) in sorted(experiments.items())],
    )
    return n_rows, len(experiments)


def _discover_datasets(expression_dir: str) -> list[str]:
    if not os.path.isdir(expression_dir):
        return []
    out: list[str] = []
    for name in sorted(os.listdir(expression_dir)):
        if name.startswith(".") or name == "README.md":
            continue
        path = os.path.join(expression_dir, name)
        if not os.path.isdir(path):
            continue
        if _find_file(path, ("meta", "samples_meta", "samples")) and _find_file(
            path, ("cpm",)
        ):
            out.append(name)
    return out


def build_expression_db(
    *,
    expression_dir: str = EXPRESSION_DIR,
    out_path: str | None = None,
    force: bool = False,
) -> str:
    datasets = _discover_datasets(expression_dir)
    if not datasets:
        raise RuntimeError(
            f"No expression datasets found under {expression_dir}/ "
            "(need <dataset>/meta.tsv + cpm.tsv)"
        )

    dest = out_path or expression_db_path()
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    if os.path.isfile(dest) and not force:
        raise RuntimeError(f"{dest} exists; pass --force to rebuild")
    if os.path.isfile(dest):
        os.remove(dest)

    con = sqlite3.connect(dest)
    try:
        con.execute("PRAGMA journal_mode=OFF")
        con.execute("PRAGMA synchronous=OFF")
        con.execute("PRAGMA temp_store=MEMORY")
        _ensure_schema(con)
        con.execute(
            "INSERT INTO meta(key, value) VALUES (?, ?)",
            ("schema_version", str(SCHEMA_VERSION)),
        )
        con.execute(
            "INSERT INTO meta(key, value) VALUES (?, ?)",
            ("built_at", time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())),
        )

        for dataset_id in datasets:
            ddir = os.path.join(expression_dir, dataset_id)
            meta_path = _find_file(ddir, ("meta", "samples_meta", "samples"))
            cpm_path = _find_file(ddir, ("cpm",))
            deg_path = _find_file(ddir, ("deg",))
            assert meta_path and cpm_path

            print(f"Dataset {dataset_id}:")
            cfg = _load_dataset_json(ddir)
            samples = _read_meta(meta_path)
            sample_accs = [r["sample_acc"] for r in samples]
            if len(sample_accs) != len(set(sample_accs)):
                raise RuntimeError(f"{meta_path}: duplicate sample_acc values")

            genome = (cfg.get("genome") or "").strip() or _majority_genome(samples)
            label = (cfg.get("label") or "").strip() or dataset_id
            if not genome:
                raise RuntimeError(
                    f"{ddir}: set genome in dataset.json or provide a genome "
                    "column in meta.tsv"
                )

            con.execute(
                "INSERT INTO datasets(dataset_id, genome, label) VALUES (?, ?, ?)",
                (dataset_id, genome, label),
            )
            con.executemany(
                """
                INSERT INTO samples(
                    dataset_id, sample_idx, sample_acc, experiment_acc, organ,
                    "group", genome, species, stage, inbred, short_title, title
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                [
                    (
                        dataset_id,
                        idx,
                        r["sample_acc"],
                        r.get("experiment_acc") or "",
                        r.get("organ") or "",
                        r.get("group") or "",
                        r.get("genome") or genome,
                        r.get("species") or "",
                        r.get("stage") or "",
                        r.get("inbred") or "",
                        r.get("short_title") or "",
                        r.get("title") or "",
                    )
                    for idx, r in enumerate(samples)
                ],
            )
            print(f"  samples: {len(samples)}  genome={genome!r}  label={label!r}")

            t0 = time.time()
            n_genes = _ingest_cpm(con, dataset_id, cpm_path, sample_accs)
            print(f"  cpm genes: {n_genes}  ({time.time() - t0:.1f}s)")

            if deg_path:
                t0 = time.time()
                n_deg, n_exp = _ingest_deg(con, dataset_id, deg_path)
                print(
                    f"  deg rows: {n_deg}  experiments: {n_exp}  ({time.time() - t0:.1f}s)"
                )
            else:
                print("  deg: (none)")

            con.commit()

        con.execute("VACUUM")
    finally:
        con.close()

    size_mb = os.path.getsize(dest) / 1e6
    print(f"Wrote {dest} ({size_mb:.1f} MB)")
    return dest


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--expression-dir",
        default=EXPRESSION_DIR,
        help="Directory of per-dataset folders (default: expression/)",
    )
    p.add_argument(
        "--out",
        default=None,
        help=f"Output SQLite path (default: expression/{EXPRESSION_DB_NAME})",
    )
    p.add_argument("--force", action="store_true", help="Overwrite existing DB")
    args = p.parse_args(argv)
    try:
        build_expression_db(
            expression_dir=args.expression_dir,
            out_path=args.out,
            force=args.force,
        )
    except RuntimeError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
