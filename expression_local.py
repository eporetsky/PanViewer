"""
Local gene-omics lookup from ``expression/expression.db``.

When this file exists, PanViewer serves tissue + DEG JSON in the same shape as
PlantApp ``GET /api/gene-omics`` (compact tissue with ``group_stats``, indexed DEG)
and does not call the remote API.
"""

from __future__ import annotations

import math
import os
import sqlite3
import struct
import zlib
from collections import defaultdict
from typing import Any

EXPRESSION_DB_NAME = "expression.db"
SCHEMA_VERSION = 1
# Stored CPM: round(cpm * CPM_SCALE) as uint16 (matches PlantApp-style 1-decimal CPM).
CPM_SCALE = 10

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
EXPRESSION_DIR = os.path.join(BASE_DIR, "expression")


def expression_db_path(base_dir: str | None = None) -> str:
    root = base_dir if base_dir is not None else BASE_DIR
    return os.path.join(root, "expression", EXPRESSION_DB_NAME)


def expression_db_available(base_dir: str | None = None) -> bool:
    path = expression_db_path(base_dir)
    return os.path.isfile(path) and os.path.getsize(path) > 0


def _connect_ro(path: str) -> sqlite3.Connection:
    con = sqlite3.connect(f"file:{os.path.abspath(path)}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    return con


def _decode_cpm_blob(blob: bytes, n_samples: int) -> list[float | None]:
    raw = zlib.decompress(blob)
    expected = n_samples * 2
    if len(raw) != expected:
        raise RuntimeError(
            f"CPM blob length {len(raw)} != expected {expected} for {n_samples} samples"
        )
    vals = struct.unpack(f"<{n_samples}H", raw)
    out: list[float | None] = []
    for v in vals:
        # 65535 reserved for missing
        if v == 65535:
            out.append(None)
        else:
            out.append(v / float(CPM_SCALE))
    return out


def _encode_cpm_values(values: list[float | None]) -> bytes:
    packed: list[int] = []
    for v in values:
        if v is None or (isinstance(v, float) and (math.isnan(v) or math.isinf(v))):
            packed.append(65535)
        else:
            packed.append(min(65534, max(0, int(round(float(v) * CPM_SCALE)))))
    return zlib.compress(struct.pack(f"<{len(packed)}H", *packed), level=6)


def encode_cpm_values(values: list[float | None]) -> bytes:
    """Public helper for ``build_expression.py``."""
    return _encode_cpm_values(values)


def _pick_dataset(
    con: sqlite3.Connection, gene_id: str, genome: str | None
) -> sqlite3.Row | None:
    """
    Resolve which dataset owns ``gene_id``.

    Prefer an exact ``datasets.genome`` match when ``genome`` is set; otherwise the
    first dataset that has the gene in ``cpm``.
    """
    gid = gene_id.strip()
    if genome:
        g = genome.strip()
        row = con.execute(
            """
            SELECT d.dataset_id, d.genome, d.label
            FROM datasets d
            JOIN cpm c ON c.dataset_id = d.dataset_id
            WHERE c.gene_id = ? AND d.genome = ?
            LIMIT 1
            """,
            (gid, g),
        ).fetchone()
        if row:
            return row
        # Also accept alias match on sample-level genome strings.
        row = con.execute(
            """
            SELECT d.dataset_id, d.genome, d.label
            FROM datasets d
            JOIN cpm c ON c.dataset_id = d.dataset_id
            WHERE c.gene_id = ?
              AND EXISTS (
                SELECT 1 FROM samples s
                WHERE s.dataset_id = d.dataset_id AND s.genome = ?
              )
            LIMIT 1
            """,
            (gid, g),
        ).fetchone()
        if row:
            return row

    return con.execute(
        """
        SELECT d.dataset_id, d.genome, d.label
        FROM datasets d
        JOIN cpm c ON c.dataset_id = d.dataset_id
        WHERE c.gene_id = ?
        ORDER BY d.dataset_id
        LIMIT 1
        """,
        (gid,),
    ).fetchone()


def _stdev(vals: list[float]) -> float:
    n = len(vals)
    if n < 2:
        return 0.0
    mean = sum(vals) / n
    var = sum((x - mean) ** 2 for x in vals) / (n - 1)
    return math.sqrt(var)


def _tissue_group_stats(
    con: sqlite3.Connection, dataset_id: str, gene_id: str
) -> dict[str, Any] | None:
    samples = con.execute(
        """
        SELECT sample_idx, sample_acc, experiment_acc, organ, "group", genome,
               species, stage, inbred, short_title, title
        FROM samples
        WHERE dataset_id = ?
        ORDER BY sample_idx
        """,
        (dataset_id,),
    ).fetchall()
    if not samples:
        return {"format": "compact", "group_stats": True, "groups": []}

    row = con.execute(
        "SELECT cpm_blob FROM cpm WHERE dataset_id = ? AND gene_id = ?",
        (dataset_id, gene_id),
    ).fetchone()
    if not row:
        return None

    cpm_vals = _decode_cpm_blob(bytes(row["cpm_blob"]), len(samples))
    buckets: dict[tuple[str, str], list[tuple[sqlite3.Row, float]]] = defaultdict(list)
    for s, v in zip(samples, cpm_vals):
        if v is None:
            continue
        key = (s["experiment_acc"] or "", s["group"] or "")
        buckets[key].append((s, float(v)))

    groups: list[dict[str, Any]] = []
    for (_ea, _grp), items in sorted(
        buckets.items(),
        key=lambda kv: (
            kv[1][0][0]["organ"] or "",
            kv[0][0],
            kv[0][1],
        ),
    ):
        meta = items[0][0]
        nums = [v for _, v in items]
        groups.append(
            {
                "experiment_acc": meta["experiment_acc"] or "",
                "genome": meta["genome"] or "",
                "group": meta["group"] or "",
                "inbred": meta["inbred"] or "",
                "mean_cpm": sum(nums) / len(nums),
                "n_samples": len(nums),
                "organ": meta["organ"] or "",
                "short_title": meta["short_title"] or "",
                "species": meta["species"] or "",
                "stage": meta["stage"] or "",
                "stdev_cpm": _stdev(nums),
                "title": meta["title"] or "",
            }
        )
    return {"format": "compact", "group_stats": True, "groups": groups}


def _deg_payload(
    con: sqlite3.Connection, dataset_id: str, gene_id: str
) -> dict[str, Any]:
    all_cats = [
        r[0]
        for r in con.execute(
            """
            SELECT DISTINCT category FROM deg
            WHERE dataset_id = ? AND category IS NOT NULL AND category != ''
            ORDER BY category
            """,
            (dataset_id,),
        ).fetchall()
    ]

    rows = con.execute(
        """
        SELECT experiment_acc, category, comparison, group1, group2, logFC, FDR
        FROM deg
        WHERE dataset_id = ? AND gene_id = ?
        ORDER BY category, experiment_acc, comparison
        """,
        (dataset_id, gene_id),
    ).fetchall()

    exp_meta = {
        r["experiment_acc"]: {
            "experiment_acc": r["experiment_acc"],
            "short_title": r["short_title"] or "",
            "title": r["title"] or "",
        }
        for r in con.execute(
            """
            SELECT experiment_acc, short_title, title
            FROM experiments
            WHERE dataset_id = ?
            """,
            (dataset_id,),
        )
    }

    exp_order: list[dict[str, str]] = []
    exp_index: dict[str, int] = {}
    records: list[dict[str, Any]] = []
    cats_with = set()

    for r in rows:
        ea = r["experiment_acc"] or ""
        if ea not in exp_index:
            exp_index[ea] = len(exp_order)
            meta = exp_meta.get(ea) or {
                "experiment_acc": ea,
                "short_title": "",
                "title": "",
            }
            exp_order.append(meta)
        cat = r["category"] or ""
        if cat:
            cats_with.add(cat)
        records.append(
            {
                "FDR": r["FDR"],
                "category": cat,
                "comparison": r["comparison"] or "",
                "e": exp_index[ea],
                "group1": r["group1"] or "",
                "group2": r["group2"] or "",
                "logFC": r["logFC"],
            }
        )

    categories = all_cats if all_cats else sorted(cats_with)
    without = [c for c in categories if c not in cats_with]

    return {
        "format": "indexed",
        "categories": categories,
        "categories_without_deg_for_gene": without,
        "experiment_order": exp_order,
        "records": records,
    }


def fetch_local_omics(
    gene_id: str, *, genome: str | None = None, base_dir: str | None = None
) -> dict[str, Any]:
    """
    Return the same dict contract as ``fetch_plantapp_omics``.
    """
    gid = (gene_id or "").strip()
    out: dict[str, Any] = {
        "query_gene_id": gid,
        "ok": False,
        "unknown_gene": False,
        "tissue": None,
        "deg": None,
        "resolved_gene_id": None,
        "resolved_genome": None,
        "source": "local",
    }
    if not gid:
        out["error"] = "gene_id is required"
        return out

    path = expression_db_path(base_dir)
    if not os.path.isfile(path):
        out["error"] = f"{EXPRESSION_DB_NAME} not found"
        return out

    try:
        con = _connect_ro(path)
    except sqlite3.Error as e:
        out["error"] = f"Cannot open {EXPRESSION_DB_NAME}: {e}"
        return out

    try:
        ds = _pick_dataset(con, gid, genome)
        if not ds:
            out["unknown_gene"] = True
            out["ok"] = True
            return out

        dataset_id = ds["dataset_id"]
        tissue = _tissue_group_stats(con, dataset_id, gid)
        deg = _deg_payload(con, dataset_id, gid)
        out["resolved_gene_id"] = gid
        out["resolved_genome"] = ds["genome"] or None
        out["tissue"] = tissue
        out["deg"] = deg
        out["ok"] = True
        return out
    except Exception as e:
        out["error"] = str(e)
        return out
    finally:
        con.close()


# Re-export for callers that only import this module
__all__ = [
    "EXPRESSION_DB_NAME",
    "SCHEMA_VERSION",
    "CPM_SCALE",
    "EXPRESSION_DIR",
    "expression_db_path",
    "expression_db_available",
    "encode_cpm_values",
    "fetch_local_omics",
]
