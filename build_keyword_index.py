"""
Build the shared keyword search index used by the /search keyword fallback.

Inputs (under ``search/``):

* ``reference/Araport11_functional_descriptions_20250331.txt`` (transcript-keyed):
  ``name``, ``gene_model_type``, ``short_description``, ``Curator_summary``,
  ``Computational_description``. Multiple transcripts per locus; we keep the
  transcript row whose ``short + curator + computational`` description is
  longest.
* ``reference/gene_aliases_20250331.txt``: ``locus_name``, ``symbol``, ``full_name``.
* ``reference/IRGSP-1.0_representative_annotation_2025-03-19.tsv``: transcript-keyed
  rice annotation. We collapse multiple transcripts per ``Locus_ID`` to the row
  with the longest ``Description``, merge curated symbol/name synonym columns
  (RAP / CGSNL / Oryzabase), and index **only** that manual text (no GO /
  InterPro columns from the file).
* ``reference/geneInfo.table.txt``: rice ``Symbol`` (``|``-separated) -> RAPdb
  (``Os01g0100100``) and MSU (``LOC_Os01g46460``); extra synonyms / aliases.
* ``mmseqs_<species>/<accession>.tsv`` (one per accession per species), e.g.
  ``mmseqs_wheat``, ``mmseqs_barley``, ``mmseqs_oat``. Legacy single-dir layout
  ``mmseqs/`` is still accepted. Required columns: ``gene_id``, ``arabidopsis``,
  ``rice``. Optional: ``arabidopsis_evalue``, ``arabidopsis_pident``,
  ``rice_evalue``, ``rice_pident`` (from MMseqs2 ``--format-output``). Best hit
  per row; each TSV is matched to the ``database/<species>.db`` that contains
  that accession.

Output: ``search/keyword_index.sqlite`` with ``refs``, ``ref_aliases``,
``refs_fts`` (FTS5), ``pan_hits``, and ``meta``.

Fail-fast (per AGENTS.md): missing required input files raise; mmseqs TSVs with
zero non-empty hit rows print a clear warning but do not silently invent rows.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import re
import sqlite3
import sys
import time
from collections import defaultdict
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parent
SEARCH_DIR = REPO_ROOT / "search"
DB_DIR = REPO_ROOT / "database"

REF_DIR_DEFAULT = SEARCH_DIR / "reference"
OUT_DEFAULT = SEARCH_DIR / "keyword_index.sqlite"

AT_DESC_FILE = "Araport11_functional_descriptions_20250331.txt"
AT_ALIASES_FILE = "gene_aliases_20250331.txt"
OS_IRGSP_FILE = "IRGSP-1.0_representative_annotation_2025-03-19.tsv"
OS_GENEINFO_FILE = "geneInfo.table.txt"

AT_LOCUS_RE = re.compile(r"^(AT[1-5MC]G\d{5})(?:\.\d+)?$", re.IGNORECASE)
OS_LOCUS_RE = re.compile(r"^(Os\d{1,2})g(\d{7})$", re.IGNORECASE)
OS_TRANSCRIPT_RE = re.compile(r"^(Os\d{1,2})t(\d{7})(?:-\d+)?$", re.IGNORECASE)
MSU_LOCUS_RE = re.compile(r"^LOC_Os\d{1,2}g\d{5}$", re.IGNORECASE)


# ---------- helpers ----------


def sha1_of(path: Path) -> str:
    h = hashlib.sha1()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def is_null(s) -> bool:
    if s is None:
        return True
    t = str(s).strip()
    if not t:
        return True
    return t.lower() in {"null", "none", "na", "n/a", "."}


def clean_field(s) -> str:
    if is_null(s):
        return ""
    return str(s).strip()


def parse_float_opt(s) -> float | None:
    """Parse scientific or decimal floats from mmseqs columns; None if empty."""
    if is_null(s):
        return None
    t = str(s).strip().replace(",", "")
    if not t:
        return None
    try:
        return float(t)
    except ValueError:
        return None


def split_synonyms(s: str) -> list[str]:
    """Split a synonym cell on common separators; drop blanks/NULLs."""
    if is_null(s):
        return []
    parts = re.split(r"\s*[,;|/]\s*", str(s))
    out: list[str] = []
    seen: set[str] = set()
    for p in parts:
        t = p.strip()
        if not t or is_null(t):
            continue
        if t in seen:
            continue
        seen.add(t)
        out.append(t)
    return out


def join_pipe(items) -> str:
    seen: set[str] = set()
    out: list[str] = []
    for it in items:
        t = (it or "").strip()
        if not t:
            continue
        if t in seen:
            continue
        seen.add(t)
        out.append(t)
    return "|".join(out)


def at_locus_from_name(name: str) -> str | None:
    s = (name or "").strip()
    if not s:
        return None
    m = AT_LOCUS_RE.match(s)
    if m:
        return m.group(1).upper()
    return None


def normalize_os_to_locus(token: str) -> str | None:
    """Accept ``Os01g0100100``, ``Os01t0100100``, ``Os01t0100100-01``. Return locus."""
    s = (token or "").strip()
    if not s:
        return None
    m = OS_LOCUS_RE.match(s)
    if m:
        return f"{m.group(1)}g{m.group(2)}"
    m = OS_TRANSCRIPT_RE.match(s)
    if m:
        return f"{m.group(1)}g{m.group(2)}"
    return None


def norm_alias(s: str) -> str:
    return (s or "").strip().lower()


# ---------- parsers ----------


def parse_at_descriptions(path: Path) -> dict[str, dict]:
    """
    Locus -> {short_description, curator_summary, computational_description,
              transcript_ids: [..]}.
    Per locus we keep the transcript row whose combined description is longest.
    """
    if not path.is_file():
        raise FileNotFoundError(f"Missing Arabidopsis descriptions file: {path}")
    rows_by_locus: dict[str, dict] = {}
    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {
            "name",
            "gene_model_type",
            "short_description",
            "Curator_summary",
            "Computational_description",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise RuntimeError(
                f"{path}: missing required columns: {sorted(missing)}"
            )
        for row in reader:
            name = (row.get("name") or "").strip()
            locus = at_locus_from_name(name)
            if not locus:
                continue
            sd = clean_field(row.get("short_description"))
            cu = clean_field(row.get("Curator_summary"))
            cp = clean_field(row.get("Computational_description"))
            combined_len = len(sd) + len(cu) + len(cp)
            cur = rows_by_locus.get(locus)
            if cur is None:
                rows_by_locus[locus] = {
                    "short_description": sd,
                    "curator_summary": cu,
                    "computational_description": cp,
                    "transcript_ids": [name],
                    "_combined_len": combined_len,
                }
            else:
                if name and name not in cur["transcript_ids"]:
                    cur["transcript_ids"].append(name)
                if combined_len > cur["_combined_len"]:
                    cur["short_description"] = sd
                    cur["curator_summary"] = cu
                    cur["computational_description"] = cp
                    cur["_combined_len"] = combined_len
    return rows_by_locus


def parse_at_aliases(path: Path) -> dict[str, dict]:
    """locus -> {symbols: set, full_names: set}."""
    if not path.is_file():
        raise FileNotFoundError(f"Missing Arabidopsis aliases file: {path}")
    out: dict[str, dict] = {}
    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"locus_name", "symbol", "full_name"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise RuntimeError(
                f"{path}: missing required columns: {sorted(missing)}"
            )
        for row in reader:
            locus = (row.get("locus_name") or "").strip().upper()
            if not locus:
                continue
            sym = clean_field(row.get("symbol"))
            fn = clean_field(row.get("full_name"))
            ent = out.setdefault(locus, {"symbols": [], "full_names": []})
            if sym and sym not in ent["symbols"]:
                ent["symbols"].append(sym)
            if fn and fn not in ent["full_names"]:
                ent["full_names"].append(fn)
    return out


def parse_irgsp(path: Path) -> dict[str, dict]:
    """
    Locus_ID -> merged annotation, picking the transcript with the longest
    ``Description``. Curated text only: symbol/name synonym columns from the
    IRGSP table (no GO / InterPro).
    """
    if not path.is_file():
        raise FileNotFoundError(f"Missing IRGSP file: {path}")
    SYM_COLS = (
        "RAP-DB Gene Symbol Synonym(s)",
        "CGSNL Gene Symbol",
        "Oryzabase Gene Symbol Synonym(s)",
    )
    NAME_COLS = (
        "RAP-DB Gene Name Synonym(s)",
        "CGSNL Gene Name",
        "Oryzabase Gene Name Synonym(s)",
    )
    out: dict[str, dict] = {}
    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        header = reader.fieldnames or []
        required = (
            "Transcript_ID",
            "Locus_ID",
            "Description",
            *SYM_COLS,
            *NAME_COLS,
        )
        missing = [c for c in required if c not in header]
        if missing:
            raise RuntimeError(
                f"{path}: missing required column(s): {missing}"
            )
        for row in reader:
            locus = (row.get("Locus_ID") or "").strip()
            if not locus:
                continue
            locus_norm = normalize_os_to_locus(locus) or locus
            tx = (row.get("Transcript_ID") or "").strip()
            desc = clean_field(row.get("Description"))
            ent = out.setdefault(
                locus_norm,
                {
                    "description": "",
                    "symbols": [],
                    "names": [],
                    "transcript_ids": [],
                    "_desc_len": -1,
                },
            )
            if tx and tx not in ent["transcript_ids"]:
                ent["transcript_ids"].append(tx)
            for col in SYM_COLS:
                for tok in split_synonyms(row.get(col, "")):
                    if tok not in ent["symbols"]:
                        ent["symbols"].append(tok)
            for col in NAME_COLS:
                for tok in split_synonyms(row.get(col, "")):
                    if tok not in ent["names"]:
                        ent["names"].append(tok)
            if len(desc) > ent["_desc_len"]:
                ent["description"] = desc
                ent["_desc_len"] = len(desc)
    return out


def parse_os_geneinfo(path: Path) -> dict[str, dict]:
    """
    RAPdb locus -> {symbols: [..], msu: [..]}. The Symbol column may carry
    pipe-separated synonyms (``ABC1|OsFd-GOGAT|SPL32|ES7``).
    """
    if not path.is_file():
        # Not strictly required, but if present it adds curated symbols + MSU.
        return {}
    out: dict[str, dict] = {}
    with open_text(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        header = reader.fieldnames or []
        for col in ("Symbol", "RAPdb", "MSU"):
            if col not in header:
                raise RuntimeError(f"{path}: missing required column: {col}")
        for row in reader:
            rapdb = (row.get("RAPdb") or "").strip()
            if not rapdb:
                continue
            locus = normalize_os_to_locus(rapdb) or rapdb
            sym_raw = (row.get("Symbol") or "").strip()
            msu = (row.get("MSU") or "").strip()
            ent = out.setdefault(locus, {"symbols": [], "msu": []})
            if sym_raw:
                for tok in [p.strip() for p in sym_raw.split("|")]:
                    if tok and not is_null(tok) and tok not in ent["symbols"]:
                        ent["symbols"].append(tok)
            if msu and msu not in ent["msu"]:
                ent["msu"].append(msu)
    return out


# ---------- DB build ----------


SCHEMA = [
    """
    CREATE TABLE refs (
        rowid       INTEGER PRIMARY KEY,
        ref_id      TEXT UNIQUE NOT NULL,
        species     TEXT NOT NULL,
        symbols     TEXT,
        full_name   TEXT,
        description TEXT,
        synonyms    TEXT,
        at_short_description      TEXT,
        at_curator_summary        TEXT,
        at_computational_description TEXT,
        source      TEXT
    )
    """,
    """
    CREATE TABLE ref_aliases (
        alias_norm TEXT NOT NULL,
        ref_id     TEXT NOT NULL,
        kind       TEXT NOT NULL,
        PRIMARY KEY (alias_norm, ref_id, kind)
    )
    """,
    "CREATE INDEX idx_ref_aliases_norm ON ref_aliases(alias_norm)",
    """
    CREATE VIRTUAL TABLE refs_fts USING fts5(
        symbols, full_name, description, synonyms,
        at_short_description, at_curator_summary, at_computational_description,
        content='refs', content_rowid='rowid',
        tokenize='unicode61 remove_diacritics 2'
    )
    """,
    """
    CREATE TABLE pan_hits (
        species    TEXT NOT NULL,
        ref_id     TEXT NOT NULL,
        dataset_id TEXT NOT NULL,
        accession  TEXT NOT NULL,
        gene_id    TEXT NOT NULL,
        pangene    TEXT,
        evalue     REAL,
        pident     REAL,
        PRIMARY KEY (dataset_id, gene_id, species, ref_id)
    )
    """,
    "CREATE INDEX idx_pan_hits_lookup  ON pan_hits(species, ref_id, dataset_id)",
    "CREATE INDEX idx_pan_hits_pangene ON pan_hits(dataset_id, pangene)",
    "CREATE TABLE meta (k TEXT PRIMARY KEY, v TEXT)",
]


def init_db(out_path: Path, *, force: bool) -> sqlite3.Connection:
    if out_path.exists():
        if not force:
            raise FileExistsError(
                f"{out_path} already exists; pass --force to rebuild."
            )
        out_path.unlink()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(out_path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    for stmt in SCHEMA:
        conn.executescript(stmt)
    return conn


def build_refs(
    conn: sqlite3.Connection,
    at_descs: dict[str, dict],
    at_aliases: dict[str, dict],
    rice_irgsp: dict[str, dict],
    rice_geneinfo: dict[str, dict],
) -> tuple[int, int, dict[str, int]]:
    refs_rows: list[tuple] = []
    alias_rows: list[tuple] = []
    seen_alias: set[tuple[str, str, str]] = set()

    def add_alias(alias: str, ref_id: str, kind: str) -> None:
        a = norm_alias(alias)
        if not a:
            return
        key = (a, ref_id, kind)
        if key in seen_alias:
            return
        seen_alias.add(key)
        alias_rows.append(key)

    at_loci = set(at_descs.keys()) | set(at_aliases.keys())
    for locus in sorted(at_loci):
        d = at_descs.get(locus, {})
        a = at_aliases.get(locus, {})
        symbols = list(a.get("symbols", []))
        full_names = list(a.get("full_names", []))
        sd = d.get("short_description", "")
        cu = d.get("curator_summary", "")
        cp = d.get("computational_description", "")
        synonyms = list(symbols) + list(full_names)
        refs_rows.append(
            (
                locus,
                "arabidopsis",
                join_pipe(symbols),
                join_pipe(full_names),
                "",
                join_pipe(synonyms),
                sd,
                cu,
                cp,
                "Araport11",
            )
        )
        add_alias(locus, locus, "id")
        for tx in d.get("transcript_ids", []) or []:
            add_alias(tx, locus, "id")
        for sym in symbols:
            add_alias(sym, locus, "symbol")
        for n in full_names:
            add_alias(n, locus, "name")

    for locus, ent in sorted(rice_irgsp.items()):
        symbols = list(ent.get("symbols", []))
        names = list(ent.get("names", []))
        gi = rice_geneinfo.get(locus, {})
        for s in gi.get("symbols", []):
            if s not in symbols:
                symbols.append(s)
        merged: list[str] = []
        for tok in symbols + names:
            if tok and tok not in merged:
                merged.append(tok)
        msu_list = list(gi.get("msu", []))
        description = ent.get("description", "")
        refs_rows.append(
            (
                locus,
                "rice",
                join_pipe(merged),
                "",
                description,
                join_pipe(msu_list),
                "",
                "",
                "",
                "IRGSP" + ("/geneInfo" if gi else ""),
            )
        )
        add_alias(locus, locus, "id")
        for tx in ent.get("transcript_ids", []) or []:
            add_alias(tx, locus, "id")
            m = OS_TRANSCRIPT_RE.match(tx)
            if m:
                add_alias(f"{m.group(1)}t{m.group(2)}", locus, "id")
        for s in merged:
            add_alias(s, locus, "symbol")
        for m in msu_list:
            add_alias(m, locus, "id")

    cur = conn.cursor()
    cur.executemany(
        "INSERT INTO refs(ref_id, species, symbols, full_name, description, "
        "synonyms, at_short_description, at_curator_summary, "
        "at_computational_description, source) VALUES (?,?,?,?,?,?,?,?,?,?)",
        refs_rows,
    )
    cur.executemany(
        "INSERT OR IGNORE INTO ref_aliases(alias_norm, ref_id, kind) VALUES (?,?,?)",
        alias_rows,
    )
    cur.execute(
        "INSERT INTO refs_fts(rowid, symbols, full_name, description, synonyms, "
        "at_short_description, at_curator_summary, at_computational_description) "
        "SELECT rowid, symbols, full_name, description, synonyms, "
        "at_short_description, at_curator_summary, at_computational_description "
        "FROM refs"
    )
    conn.commit()

    by_species: dict[str, int] = defaultdict(int)
    for r in refs_rows:
        by_species[r[1]] += 1
    return len(refs_rows), len(alias_rows), dict(by_species)


# ---------- mmseqs -> pan_hits ----------


def discover_datasets(db_dir: Path) -> list[tuple[str, Path]]:
    if not db_dir.is_dir():
        return []
    out: list[tuple[str, Path]] = []
    for p in sorted(db_dir.iterdir()):
        if p.suffix == ".db" and not p.name.startswith("_") and not p.name.startswith("."):
            out.append((p.stem, p))
    return out


def find_dataset_for_accession(
    datasets: list[tuple[str, Path]], accession: str
) -> tuple[str, Path] | None:
    """Match mmseqs TSV stem to a species DB via ``genes`` or ``gene_coords``."""
    for dsid, dbp in datasets:
        try:
            with sqlite3.connect(f"file:{dbp}?mode=ro", uri=True) as conn:
                cur = conn.cursor()
                for table in ("genes", "gene_coords"):
                    try:
                        cur.execute(
                            f"SELECT 1 FROM {table} WHERE accession = ? LIMIT 1",
                            (accession,),
                        )
                        if cur.fetchone() is not None:
                            return (dsid, dbp)
                    except sqlite3.Error:
                        continue
        except sqlite3.Error:
            continue
    return None


def load_gene_to_pangene(db_path: Path, accession: str) -> dict[str, str]:
    out: dict[str, str] = {}
    with sqlite3.connect(f"file:{db_path}?mode=ro", uri=True) as conn:
        conn.row_factory = sqlite3.Row
        cur = conn.cursor()
        cur.execute(
            "SELECT gene_id, pangene FROM genes WHERE accession = ?",
            (accession,),
        )
        for r in cur.fetchall():
            out[r["gene_id"]] = r["pangene"] or ""
    return out


def ingest_mmseqs(
    conn: sqlite3.Connection,
    mmseqs_dir: Path,
    datasets: list[tuple[str, Path]],
) -> dict:
    if not mmseqs_dir.is_dir():
        raise FileNotFoundError(f"Missing mmseqs directory: {mmseqs_dir}")
    cur = conn.cursor()
    stats: dict[str, dict] = {}
    files = sorted(p for p in mmseqs_dir.glob("*.tsv") if p.is_file())
    if not files:
        print(f"  No mmseqs TSVs in {mmseqs_dir}", file=sys.stderr)
        return {"files": 0, "hits": 0}
    total_hits = 0
    for tsv in files:
        accession = tsv.stem
        match = find_dataset_for_accession(datasets, accession)
        if match is None:
            print(
                f"  Skip {tsv.name}: no database/*.db has accession='{accession}'",
                file=sys.stderr,
            )
            stats[accession] = {"dataset": None, "rows": 0, "hits": 0, "missing_in_db": 0}
            continue
        dataset_id, db_path = match
        gene_to_pan = load_gene_to_pangene(db_path, accession)
        rows_in = 0
        hits_in = 0
        missing_in_db = 0
        seen_at_empty = True
        seen_os_empty = True
        batch: list[tuple] = []
        with tsv.open("r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            if "gene_id" not in (reader.fieldnames or []):
                raise RuntimeError(f"{tsv}: missing 'gene_id' column")
            for row in reader:
                rows_in += 1
                gid = (row.get("gene_id") or "").strip()
                if not gid:
                    continue
                at = (row.get("arabidopsis") or "").strip()
                os_ = (row.get("rice") or "").strip()
                if at:
                    seen_at_empty = False
                if os_:
                    seen_os_empty = False
                if not at and not os_:
                    continue
                pangene = gene_to_pan.get(gid)
                if pangene is None:
                    missing_in_db += 1
                    continue
                at_ev = parse_float_opt(row.get("arabidopsis_evalue"))
                at_pi = parse_float_opt(row.get("arabidopsis_pident"))
                os_ev = parse_float_opt(row.get("rice_evalue"))
                os_pi = parse_float_opt(row.get("rice_pident"))
                if at:
                    locus = at_locus_from_name(at)
                    if locus:
                        batch.append((
                            "arabidopsis",
                            locus,
                            dataset_id,
                            accession,
                            gid,
                            pangene,
                            at_ev,
                            at_pi,
                        ))
                        hits_in += 1
                if os_:
                    locus_os = normalize_os_to_locus(os_)
                    if locus_os:
                        batch.append((
                            "rice",
                            locus_os,
                            dataset_id,
                            accession,
                            gid,
                            pangene,
                            os_ev,
                            os_pi,
                        ))
                        hits_in += 1
                if len(batch) >= 5000:
                    cur.executemany(
                        "INSERT OR REPLACE INTO pan_hits(species, ref_id, dataset_id, "
                        "accession, gene_id, pangene, evalue, pident) "
                        "VALUES (?,?,?,?,?,?,?,?)",
                        batch,
                    )
                    batch.clear()
        if batch:
            cur.executemany(
                "INSERT OR REPLACE INTO pan_hits(species, ref_id, dataset_id, "
                "accession, gene_id, pangene, evalue, pident) "
                "VALUES (?,?,?,?,?,?,?,?)",
                batch,
            )
        conn.commit()
        total_hits += hits_in
        stats[accession] = {
            "dataset": dataset_id,
            "rows": rows_in,
            "hits": hits_in,
            "missing_in_db": missing_in_db,
        }
        warn = ""
        if seen_at_empty and seen_os_empty:
            warn = " WARN: file has no non-empty AT/OS hits"
        elif missing_in_db:
            warn = f" ({missing_in_db} gene_ids not in {dataset_id}.db)"
        print(
            f"  {tsv.name}: {hits_in} hits -> {dataset_id} (rows={rows_in}){warn}"
        )
    return {"files": len(files), "hits": total_hits, "per_file": stats}


def discover_mmseqs_dirs(search_dir: Path) -> list[Path]:
    """
    Default mmseqs inputs: every ``search/mmseqs_*`` directory, else legacy
    ``search/mmseqs`` if it exists.
    """
    dirs = sorted(
        p
        for p in search_dir.iterdir()
        if p.is_dir() and p.name.startswith("mmseqs_")
    )
    if dirs:
        return dirs
    legacy = search_dir / "mmseqs"
    return [legacy] if legacy.is_dir() else []


def ingest_mmseqs_all(
    conn: sqlite3.Connection,
    mmseqs_dirs: list[Path],
    datasets: list[tuple[str, Path]],
) -> dict:
    combined: dict = {"files": 0, "hits": 0, "dirs": []}
    for mm_dir in mmseqs_dirs:
        print(f"\nmmseqs: {mm_dir}")
        st = ingest_mmseqs(conn, mm_dir, datasets)
        combined["files"] += st.get("files", 0)
        combined["hits"] += st.get("hits", 0)
        combined["dirs"].append({"path": str(mm_dir), **st})
    return combined


# ---------- main ----------


def write_meta(conn: sqlite3.Connection, meta: dict) -> None:
    cur = conn.cursor()
    cur.executemany(
        "INSERT OR REPLACE INTO meta(k, v) VALUES (?, ?)",
        [(k, str(v)) for k, v in meta.items()],
    )
    conn.commit()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ref-dir", type=Path, default=REF_DIR_DEFAULT)
    ap.add_argument(
        "--mmseqs-dir",
        type=Path,
        action="append",
        dest="mmseqs_dirs",
        metavar="DIR",
        help=(
            "Per-accession mmseqs TSV directory (repeat for multiple species). "
            "Default: all search/mmseqs_* dirs, else search/mmseqs/"
        ),
    )
    ap.add_argument("--db-dir", type=Path, default=DB_DIR)
    ap.add_argument("--out", type=Path, default=OUT_DEFAULT)
    ap.add_argument("--force", action="store_true", help="Overwrite existing output")
    args = ap.parse_args()

    t0 = time.time()
    at_desc_path = args.ref_dir / AT_DESC_FILE
    at_alias_path = args.ref_dir / AT_ALIASES_FILE
    os_irgsp_path = args.ref_dir / OS_IRGSP_FILE
    os_geneinfo_path = args.ref_dir / OS_GENEINFO_FILE

    print(f"AT descriptions: {at_desc_path}")
    at_descs = parse_at_descriptions(at_desc_path)
    print(f"  loci: {len(at_descs)}")

    print(f"AT aliases: {at_alias_path}")
    at_aliases = parse_at_aliases(at_alias_path)
    print(f"  loci: {len(at_aliases)}")

    print(f"OS IRGSP: {os_irgsp_path}")
    rice_irgsp = parse_irgsp(os_irgsp_path)
    print(f"  loci: {len(rice_irgsp)}")

    print(f"OS geneInfo: {os_geneinfo_path}")
    rice_gi = parse_os_geneinfo(os_geneinfo_path)
    print(f"  loci with curated symbols: {len(rice_gi)}")

    conn = init_db(args.out, force=args.force)
    try:
        n_refs, n_alias, by_species = build_refs(
            conn, at_descs, at_aliases, rice_irgsp, rice_gi
        )
        print(
            f"refs={n_refs} ({by_species}), aliases={n_alias}"
        )

        datasets = discover_datasets(args.db_dir)
        print(f"datasets: {[d for d, _ in datasets]}")
        mmseqs_dirs = args.mmseqs_dirs or discover_mmseqs_dirs(SEARCH_DIR)
        if not mmseqs_dirs:
            raise FileNotFoundError(
                "No mmseqs directories found under search/. "
                "Expected mmseqs_wheat/, mmseqs_barley/, mmseqs_oat/, or legacy mmseqs/."
            )
        mm_stats = ingest_mmseqs_all(conn, mmseqs_dirs, datasets)

        meta = {
            "schema_version": 3,
            "built_at": int(time.time()),
            "at_desc_sha1": sha1_of(at_desc_path) if at_desc_path.is_file() else "",
            "at_alias_sha1": sha1_of(at_alias_path) if at_alias_path.is_file() else "",
            "os_irgsp_sha1": sha1_of(os_irgsp_path) if os_irgsp_path.is_file() else "",
            "os_geneinfo_sha1": (
                sha1_of(os_geneinfo_path) if os_geneinfo_path.is_file() else ""
            ),
            "refs_count": n_refs,
            "aliases_count": n_alias,
            "pan_hits_files": mm_stats.get("files", 0),
            "pan_hits_total": mm_stats.get("hits", 0),
            "mmseqs_dirs": ",".join(str(p) for p in mmseqs_dirs),
        }
        write_meta(conn, meta)
        conn.execute("ANALYZE")
        conn.commit()
        # Merge WAL so ``schema_version`` in meta is visible without -wal (e.g. Docker COPY).
        conn.execute("PRAGMA wal_checkpoint(TRUNCATE)")
    finally:
        conn.close()

    dt = time.time() - t0
    print(f"\nWrote {args.out} ({dt:.1f}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
