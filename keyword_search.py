"""
Runtime helpers for the keyword-search fallback in ``/search``.

The shared index at ``search/keyword_index.sqlite`` is built by
``build_keyword_index.py``. This module reads it (read-only) and returns
results shaped to drop straight into the existing ``results.html`` flow.

Pan-gene ordering (when ``pan_hits`` has ``evalue`` / ``pident`` columns):
best (minimum) e-value across retained hits, then hit density (fraction of
cluster genes with at least one retained hit). Rows with NULL ``evalue`` are
kept for backward compatibility and ignored by the optional e-value cutoff.

Fail-fast (per AGENTS.md): if the index file is missing or unreadable, the
caller can detect that via :func:`keyword_index_available` and surface a
clear error rather than silently returning empty results.
"""
from __future__ import annotations

import os
import re
import sqlite3
from typing import Any


# Tokens we send into FTS5: letters/digits, plus simple "-" inside identifiers
# kept by replacing it with a space (FTS5 treats "-" as a unary NOT operator
# unless quoted, and a too-clever escape only invites surprises).
_TOKEN_RE = re.compile(r"[A-Za-z0-9]{2,}")
# Max FTS hits and final ref shortlist size. Kept conservative so the UI stays
# fast for very generic queries like "kinase".
DEFAULT_REF_LIMIT = 25
DEFAULT_PANGENE_LIMIT = 100
# Drop mmseqs hits weaker than this e-value when the column is present (None =
# keep all rows). Rows with NULL evalue are never dropped.
DEFAULT_MAX_HIT_EVALUE = 1e-5
# Must match ``schema_version`` written by ``build_keyword_index.py``.
KEYWORD_INDEX_SCHEMA_VERSION = 3


def keyword_index_path() -> str:
    """Default index path; matches ``build_keyword_index.py`` --out default."""
    base = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(base, "search", "keyword_index.sqlite")


def keyword_index_available() -> bool:
    p = keyword_index_path()
    return os.path.isfile(p) and os.path.getsize(p) > 0


def _connect_ro() -> sqlite3.Connection:
    p = keyword_index_path()
    if not os.path.isfile(p):
        raise RuntimeError(
            f"Keyword search index not found at {p}. "
            "Run `python build_keyword_index.py` to create it."
        )
    conn = sqlite3.connect(f"file:{p}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    # If the index was built with WAL, meta may live in -wal until checkpointed.
    try:
        conn.execute("PRAGMA wal_checkpoint(PASSIVE)")
    except sqlite3.Error:
        pass
    return conn


def _keyword_index_schema_version(conn: sqlite3.Connection) -> int:
    cur = conn.cursor()
    cur.execute("SELECT v FROM meta WHERE k = 'schema_version'")
    row = cur.fetchone()
    if not row or row[0] is None:
        return 0
    try:
        return int(str(row[0]).strip())
    except ValueError:
        return 0


def _tokenize(query: str) -> list[str]:
    return [m.group(0) for m in _TOKEN_RE.finditer(query or "")]


def _build_fts_match(tokens: list[str]) -> str | None:
    """
    AND across tokens, prefix on each. ``kinase resistance`` becomes
    ``"kinase"* "resistance"*``. Tokens are already alnum-only, so quoting
    them as bare phrases is safe.
    """
    parts: list[str] = []
    for t in tokens:
        if len(t) < 2:
            continue
        parts.append(f'"{t}"*')
    if not parts:
        return None
    return " ".join(parts)


def _alias_terms(query: str, tokens: list[str]) -> list[str]:
    """Lowercased phrases to try as exact alias matches: full query + each token."""
    out: list[str] = []
    seen: set[str] = set()
    full = (query or "").strip().lower()
    if full and full not in seen:
        seen.add(full)
        out.append(full)
    for t in tokens:
        a = t.lower()
        if a and a not in seen:
            seen.add(a)
            out.append(a)
    return out


def _join_description_parts_dot(*parts: str) -> str:
    """
    Join non-empty description fragments in order; insert ``. `` between
    fragments when the accumulated text does not already end with ``.``.
    """
    out = ""
    for raw in parts:
        t = (raw or "").strip()
        if not t:
            continue
        if not out:
            out = t
            continue
        sep = " " if out.rstrip().endswith(".") else ". "
        out = f"{out}{sep}{t}"
    return out


def _snippet(desc: str, tokens: list[str], width: int = 140) -> str:
    """Lightweight snippet around the first matched token; no markup added."""
    if not desc:
        return ""
    if not tokens:
        return desc[:width] + ("..." if len(desc) > width else "")
    lower = desc.lower()
    for t in tokens:
        i = lower.find(t.lower())
        if i >= 0:
            start = max(0, i - width // 3)
            end = min(len(desc), start + width)
            prefix = "..." if start > 0 else ""
            suffix = "..." if end < len(desc) else ""
            return f"{prefix}{desc[start:end]}{suffix}"
    return desc[:width] + ("..." if len(desc) > width else "")


def _pan_hits_has_score_columns(cur: sqlite3.Cursor) -> bool:
    cur.execute("PRAGMA table_info(pan_hits)")
    names = {r[1] for r in cur.fetchall()}
    return "evalue" in names and "pident" in names


def _fmt_evalue(x: float | None) -> str:
    if x is None:
        return "—"
    if x == 0.0:
        return "0"
    ax = abs(x)
    if ax < 1e-300:
        return f"{x:.1e}"
    if ax >= 0.01:
        return f"{x:.4g}"
    return f"{x:.2e}"


def _fmt_pident(x: float | None) -> str:
    if x is None:
        return "—"
    return f"{x:.1f}%"


def _pangenes_for_gene_hits(
    db_path: str | None, rows: list[dict[str, Any]]
) -> dict[tuple[str, str], str]:
    """
    Map ``(gene_id, accession)`` → ``pangene`` in the active variant DB.

    Used when keyword ``pan_hits`` were built for another index (e.g. shared
    ``dataset_id=oat`` mmseqs) but results should be grouped by the selected
    variant (e.g. ``panoat.db``).
    """
    if not db_path or not os.path.isfile(db_path) or not rows:
        return {}
    keys = list(dict.fromkeys((str(r["gene_id"]), str(r["accession"])) for r in rows))
    out: dict[tuple[str, str], str] = {}
    chunk = 400
    uri = f"file:{os.path.abspath(db_path)}?mode=ro"
    with sqlite3.connect(uri, uri=True) as db:
        for i in range(0, len(keys), chunk):
            part = keys[i : i + chunk]
            if not part:
                continue
            ph_pairs = ",".join("(?,?)" for _ in part)
            params: list[str] = []
            for gid, acc in part:
                params.extend([gid, acc])
            cur = db.execute(
                f"SELECT gene_id, accession, pangene FROM genes "
                f"WHERE (gene_id, accession) IN ({ph_pairs})",
                params,
            )
            for r in cur.fetchall():
                pan = (r[2] or "").strip()
                if pan:
                    out[(str(r[0]), str(r[1]))] = pan
    return out


def _load_pangene_sizes(db_path: str, pangenes: list[str]) -> dict[str, int]:
    """Member-gene counts per pan-gene from the species SQLite (``genes``)."""
    out: dict[str, int] = {}
    if not db_path or not os.path.isfile(db_path):
        return out
    chunk = 400
    uri = f"file:{db_path}?mode=ro"
    with sqlite3.connect(uri, uri=True) as db:
        for i in range(0, len(pangenes), chunk):
            part = pangenes[i : i + chunk]
            if not part:
                continue
            ph = ",".join("?" * len(part))
            cur = db.execute(
                f"SELECT pangene, COUNT(*) AS n FROM genes WHERE pangene IN ({ph}) "
                "GROUP BY pangene",
                part,
            )
            for r in cur:
                out[str(r[0])] = int(r[1])
    return out


def search_keywords(
    query: str,
    dataset_id: str,
    *,
    db_path: str | None = None,
    ref_limit: int = DEFAULT_REF_LIMIT,
    pangene_limit: int = DEFAULT_PANGENE_LIMIT,
    max_hit_evalue: float | None = DEFAULT_MAX_HIT_EVALUE,
    restrict_refs: set[str] | None = None,
) -> dict[str, Any]:
    """
    Run a keyword search against the shared index.

    Returns ``{"refs": [...], "pangene_groups": [...]}`` where ``pangene_groups``
    matches the gene-ID search shape, plus per-group ``keyword_stats`` when
    score columns exist in ``pan_hits``.

    ``db_path``: SQLite path for the active dataset (``database/<id>.db``).
    Used to compute cluster size and hit density. If missing, density falls
    back to ``matched_genes / max(matched_genes, 1)``.

    ``max_hit_evalue``: if not ``None``, drop hits whose ``evalue`` is strictly
    greater than this threshold (hits with NULL ``evalue`` are kept).

    ``restrict_refs``: optional case-insensitive set of ``ref_id`` values. When
    set, only ``pan_hits`` rows whose ``ref_id`` matches one of these are used
    to build ``pangene_groups``. The full matched ``refs`` list is still
    returned (each entry carries a ``selected`` flag) so the UI can offer the
    user a chance to refine the selection without losing context.
    """
    q = (query or "").strip()
    if not q:
        return {"refs": [], "pangene_groups": []}
    ds = (dataset_id or "").strip().lower()
    if not ds:
        return {"refs": [], "pangene_groups": []}

    tokens = _tokenize(q)
    fts_expr = _build_fts_match(tokens)
    aliases = _alias_terms(q, tokens)
    restrict_norm: set[str] | None = None
    if restrict_refs is not None:
        restrict_norm = {str(x).strip().upper() for x in restrict_refs if str(x).strip()}
        if not restrict_norm:
            restrict_norm = None

    conn = _connect_ro()
    try:
        cur = conn.cursor()
        idx_ver = _keyword_index_schema_version(conn)
        if idx_ver < KEYWORD_INDEX_SCHEMA_VERSION:
            raise RuntimeError(
                f"Keyword index at {keyword_index_path()} is schema version "
                f"{idx_ver}; this app expects {KEYWORD_INDEX_SCHEMA_VERSION}. "
                "Rebuild with: python3 build_keyword_index.py --force"
            )
        has_scores = _pan_hits_has_score_columns(cur)

        ref_order: list[str] = []
        ref_seen: set[str] = set()
        ref_match_kind: dict[str, str] = {}

        if aliases:
            placeholders = ",".join("?" for _ in aliases)
            cur.execute(
                f"SELECT DISTINCT ref_id FROM ref_aliases "
                f"WHERE alias_norm IN ({placeholders})",
                aliases,
            )
            for r in cur.fetchall():
                rid = r["ref_id"]
                if rid not in ref_seen:
                    ref_seen.add(rid)
                    ref_order.append(rid)
                    ref_match_kind[rid] = "alias"

        if fts_expr:
            fts_cap = max(ref_limit * 4, 50)
            try:
                cur.execute(
                    """
                    SELECT r.ref_id AS ref_id, bm25(refs_fts) AS rank
                    FROM refs_fts JOIN refs r ON r.rowid = refs_fts.rowid
                    WHERE refs_fts MATCH ?
                    ORDER BY rank LIMIT ?
                    """,
                    (fts_expr, fts_cap),
                )
                for r in cur.fetchall():
                    rid = r["ref_id"]
                    if rid not in ref_seen:
                        ref_seen.add(rid)
                        ref_order.append(rid)
                        ref_match_kind[rid] = "fts"
            except sqlite3.OperationalError:
                pass

        if not ref_order:
            return {"refs": [], "pangene_groups": []}

        placeholders = ",".join("?" for _ in ref_order)
        cur.execute(
            f"""
            SELECT DISTINCT ref_id FROM pan_hits
            WHERE dataset_id = ? AND ref_id IN ({placeholders})
            """,
            (ds, *ref_order),
        )
        actionable: set[str] = {r["ref_id"] for r in cur.fetchall()}
        actionable_order = [rid for rid in ref_order if rid in actionable]
        actionable_order = actionable_order[:ref_limit]

        if restrict_norm is not None:
            norm_list = list(restrict_norm)
            placeholders2 = ",".join("?" for _ in norm_list)
            cur.execute(
                f"SELECT DISTINCT ref_id FROM pan_hits "
                f"WHERE dataset_id = ? AND UPPER(ref_id) IN ({placeholders2})",
                (ds, *norm_list),
            )
            for r in cur.fetchall():
                rid = r["ref_id"]
                if rid not in ref_seen:
                    ref_seen.add(rid)
                    ref_order.append(rid)
                    ref_match_kind[rid] = "selected"
                if rid not in actionable_order:
                    actionable_order.append(rid)

        if not actionable_order:
            return {"refs": [], "pangene_groups": []}

        ph = ",".join("?" for _ in actionable_order)
        cur.execute(
            f"SELECT ref_id, species, symbols, full_name, description, synonyms, "
            f"at_short_description, at_curator_summary, at_computational_description "
            f"FROM refs WHERE ref_id IN ({ph})",
            actionable_order,
        )
        ref_meta: dict[str, sqlite3.Row] = {r["ref_id"]: r for r in cur.fetchall()}

        refs_out: list[dict[str, Any]] = []
        for rid in actionable_order:
            m = ref_meta.get(rid)
            if not m:
                continue
            sym_raw = m["symbols"] or ""
            sym_list = [s.strip() for s in sym_raw.split("|") if s and s.strip()]
            fn_raw = (m["full_name"] or "").strip()
            if fn_raw:
                for s in fn_raw.split("|"):
                    t = s.strip()
                    if t and t not in sym_list:
                        sym_list.append(t)
            at_s = (m["at_short_description"] or "").strip()
            at_c = (m["at_curator_summary"] or "").strip()
            at_p = (m["at_computational_description"] or "").strip()
            desc_plain = (m["description"] or "").strip()
            if (m["species"] or "").strip().lower() == "arabidopsis":
                combined_description = _join_description_parts_dot(at_s, at_c, at_p)
                desc_for_snippet = combined_description
            else:
                combined_description = desc_plain
                desc_for_snippet = desc_plain
            refs_out.append(
                {
                    "ref_id": rid,
                    "species": m["species"],
                    "symbols": sym_raw,
                    "symbols_list": sym_list,
                    "full_name": (m["full_name"] or "").strip(),
                    "description": desc_plain,
                    "combined_description": combined_description,
                    "at_short_description": at_s,
                    "at_curator_summary": at_c,
                    "at_computational_description": at_p,
                    "snippet": _snippet(desc_for_snippet, tokens),
                    "matched_kind": ref_match_kind.get(rid, "fts"),
                    "selected": (
                        restrict_norm is not None
                        and rid.upper() in restrict_norm
                    ),
                }
            )

        # When the user has not selected any refs (no ?refs= param), skip the
        # potentially expensive pangene rollup: the UI shows only the refs
        # picker until the user submits a selection.
        if restrict_norm is None:
            return {"refs": refs_out, "pangene_groups": []}

        ref_filter_keys: set[str] | None = {
            rid for rid in actionable_order if rid.upper() in restrict_norm
        }
        if not ref_filter_keys:
            return {"refs": refs_out, "pangene_groups": []}

        if has_scores:
            cur.execute(
                f"""
                SELECT ref_id, accession, gene_id, pangene, evalue, pident
                FROM pan_hits
                WHERE dataset_id = ? AND ref_id IN ({ph})
                """,
                (ds, *actionable_order),
            )
        else:
            cur.execute(
                f"""
                SELECT ref_id, accession, gene_id, pangene,
                       NULL AS evalue, NULL AS pident
                FROM pan_hits
                WHERE dataset_id = ? AND ref_id IN ({ph})
                """,
                (ds, *actionable_order),
            )

        raw_rows: list[dict[str, Any]] = []
        for r in cur.fetchall():
            if ref_filter_keys is not None and r["ref_id"] not in ref_filter_keys:
                continue
            ev = r["evalue"]
            pid = r["pident"]
            if (
                max_hit_evalue is not None
                and has_scores
                and ev is not None
                and float(ev) > max_hit_evalue
            ):
                continue
            raw_rows.append(
                {
                    "ref_id": r["ref_id"],
                    "accession": r["accession"],
                    "gene_id": r["gene_id"],
                    "pangene": r["pangene"] or "",
                    "evalue": float(ev) if ev is not None else None,
                    "pident": float(pid) if pid is not None else None,
                }
            )

        use_db_pangene = bool(db_path and os.path.isfile(db_path))
        pan_by_gene = (
            _pangenes_for_gene_hits(db_path, raw_rows) if use_db_pangene else {}
        )

        pangene_groups: dict[str, dict[str, Any]] = {}
        pan_order: list[str] = []
        for row in raw_rows:
            gid = row["gene_id"]
            acc = row["accession"]
            if use_db_pangene:
                pan = pan_by_gene.get((gid, acc), "")
                if not pan:
                    continue
            else:
                pan = row["pangene"]
            if not pan:
                continue
            g = pangene_groups.get(pan)
            if g is None:
                g = {
                    "pangene": pan,
                    "genes": [],
                    "_gene_keys": set(),
                    "matched_refs": [],
                    "_ref_set": set(),
                    "_rows": [],
                }
                pangene_groups[pan] = g
                pan_order.append(pan)
            key = (gid, acc)
            if key not in g["_gene_keys"]:
                g["_gene_keys"].add(key)
                g["genes"].append({"gene_id": gid, "accession": acc})
            rid = row["ref_id"]
            if rid not in g["_ref_set"]:
                g["_ref_set"].add(rid)
                g["matched_refs"].append(rid)
            g["_rows"].append(row)

        sizes = _load_pangene_sizes(db_path or "", pan_order)

        def sort_key(pan: str) -> tuple:
            g = pangene_groups[pan]
            rows = g["_rows"]
            genes_matched = len(g["_gene_keys"])
            known_sz = sizes.get(pan)
            denom = float(known_sz) if known_sz and known_sz > 0 else float(
                max(genes_matched, 1)
            )
            density = genes_matched / denom
            evals = [x["evalue"] for x in rows if x["evalue"] is not None]
            min_ev = min(evals) if evals else None
            min_for_sort = min_ev if min_ev is not None else 1.0
            return (min_for_sort, -density, -genes_matched, pan)

        pan_order.sort(key=sort_key)
        pan_order = pan_order[:pangene_limit]

        out_groups: list[dict[str, Any]] = []
        for pan in pan_order:
            g = pangene_groups[pan]
            rows = g["_rows"]
            g.pop("_gene_keys", None)
            g.pop("_ref_set", None)
            g.pop("_rows", None)

            acc_hits: dict[str, int] = {}
            for x in rows:
                acc_hits[x["accession"]] = acc_hits.get(x["accession"], 0) + 1
            hits_by_genome = sorted(
                [{"accession": a, "hits": n} for a, n in acc_hits.items() if n > 0],
                key=lambda t: (-t["hits"], t["accession"]),
            )
            genomes_with_hits = len(hits_by_genome)
            total_rows = sum(t["hits"] for t in hits_by_genome)

            genes_matched = len(g["genes"])
            known_sz = sizes.get(pan)
            denom = float(known_sz) if known_sz and known_sz > 0 else float(
                max(genes_matched, 1)
            )
            density_pct = round(100.0 * genes_matched / denom, 2)

            evals = [x["evalue"] for x in rows if x["evalue"] is not None]
            pidents = [x["pident"] for x in rows if x["pident"] is not None]
            min_ev = min(evals) if evals else None
            max_pi = max(pidents) if pidents else None

            g["keyword_stats"] = {
                "pangene_size": known_sz,
                "matched_gene_count": genes_matched,
                "hit_density_pct": density_pct,
                "genomes_with_hits": genomes_with_hits,
                "total_hit_rows": total_rows,
                "hits_by_genome": hits_by_genome,
                "min_evalue": min_ev,
                "max_pident": max_pi,
                "min_evalue_label": _fmt_evalue(min_ev),
                "max_pident_label": _fmt_pident(max_pi),
                "has_scores": has_scores,
            }
            out_groups.append(g)

        return {"refs": refs_out, "pangene_groups": out_groups}
    finally:
        conn.close()


def _annotation_text_for_ref_row(species: str, m: sqlite3.Row) -> str:
    sp = (species or "").strip().lower()
    if sp == "arabidopsis":
        at_s = (m["at_short_description"] or "").strip()
        at_c = (m["at_curator_summary"] or "").strip()
        at_p = (m["at_computational_description"] or "").strip()
        combined = _join_description_parts_dot(at_s, at_c, at_p)
        return combined or (m["description"] or "").strip()
    return (m["description"] or "").strip()


def _ev_sort_tuple(ev: Any) -> tuple[int, float]:
    """Lower tuple is better e-value; NULL sorts last."""
    if ev is None:
        return (1, 0.0)
    try:
        return (0, float(ev))
    except (TypeError, ValueError):
        return (1, 0.0)


def _prefer_hit_row(a: sqlite3.Row, b: sqlite3.Row) -> sqlite3.Row:
    """Pick the better of two pan_hits+refs rows for the same gene and species."""
    ta, va = _ev_sort_tuple(a["evalue"])
    tb, vb = _ev_sort_tuple(b["evalue"])
    if ta != tb:
        return a if ta < tb else b
    if va != vb:
        return a if va < vb else b
    pa = a["pident"]
    pb = b["pident"]
    if pa is not None and pb is not None and pa != pb:
        return a if float(pa) > float(pb) else b
    if pa is not None and pb is None:
        return a
    if pb is not None and pa is None:
        return b
    return a if str(a["ref_id"]) <= str(b["ref_id"]) else b


def cross_species_best_annotations_for_genes(
    dataset_id: str,
    gene_ids: list[str],
) -> dict[str, dict[str, Any | None]]:
    """
    Best Arabidopsis and rice keyword-index annotation per member gene (by
    minimum ``evalue`` in ``pan_hits``, tie-break higher ``pident``).

    Returns ``{gene_id: {"arabidopsis": None | {...}, "rice": None | {...}}}``
    with keys ``ref_id``, ``symbols``, ``annotation``, ``evalue``, ``pident``,
    ``evalue_label``, ``pident_label``.
    """
    ds = (dataset_id or "").strip().lower()
    uniq: list[str] = list(dict.fromkeys(str(g).strip() for g in gene_ids if str(g).strip()))
    if not ds or not uniq:
        return {}
    if not keyword_index_available():
        return {}
    out: dict[str, dict[str, Any | None]] = {}
    chunk = 400
    conn = _connect_ro()
    try:
        if _keyword_index_schema_version(conn) < KEYWORD_INDEX_SCHEMA_VERSION:
            return {}
        cur = conn.cursor()
        for i in range(0, len(uniq), chunk):
            part = uniq[i : i + chunk]
            ph = ",".join("?" * len(part))
            cur.execute(
                f"""
                SELECT h.gene_id AS gene_id, h.species AS species, h.ref_id AS ref_id,
                       h.evalue AS evalue, h.pident AS pident,
                       r.symbols AS symbols, r.full_name AS full_name,
                       r.description AS description,
                       r.at_short_description AS at_short_description,
                       r.at_curator_summary AS at_curator_summary,
                       r.at_computational_description AS at_computational_description
                FROM pan_hits h
                JOIN refs r ON r.ref_id = h.ref_id
                WHERE h.dataset_id = ? AND h.gene_id IN ({ph})
                  AND lower(h.species) IN ('arabidopsis', 'rice')
                """,
                (ds, *part),
            )
            best: dict[tuple[str, str], sqlite3.Row] = {}
            for row in cur.fetchall():
                gid = str(row["gene_id"])
                sp = str(row["species"] or "").strip().lower()
                if sp not in ("arabidopsis", "rice"):
                    continue
                key = (gid, sp)
                prev = best.get(key)
                best[key] = row if prev is None else _prefer_hit_row(prev, row)
            for (gid, sp), row in best.items():
                if gid not in out:
                    out[gid] = {"arabidopsis": None, "rice": None}
                ann = _annotation_text_for_ref_row(sp, row)
                ev = row["evalue"]
                pid = row["pident"]
                out[gid][sp] = {
                    "ref_id": row["ref_id"],
                    "symbols": (row["symbols"] or "").strip(),
                    "annotation": ann,
                    "evalue": float(ev) if ev is not None else None,
                    "pident": float(pid) if pid is not None else None,
                    "evalue_label": _fmt_evalue(float(ev) if ev is not None else None),
                    "pident_label": _fmt_pident(float(pid) if pid is not None else None),
                }
    finally:
        conn.close()
    for gid in uniq:
        out.setdefault(gid, {"arabidopsis": None, "rice": None})
    return out
