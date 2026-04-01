import csv
import hashlib
import io
import json
import logging
import math
import os
import re
import shutil
import sqlite3
import subprocess
import tempfile
import threading
import zlib
from collections import defaultdict
from typing import Any

from markupsafe import escape

from plantapp_omics import fetch_plantapp_omics

from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from flask import (
    Flask,
    Response,
    jsonify,
    render_template,
    request,
)

app = Flask(__name__)

# Always serve under /panviewer so links work at https://graingenes.org/panviewer/
# and locally at http://localhost:5050/panviewer/
APPLICATION_ROOT = "/panviewer"

class PrefixMiddleware:
    def __init__(self, app, prefix):
        self.app = app
        self.prefix = prefix.rstrip("/") if prefix else ""

    def __call__(self, environ, start_response):
        if not self.prefix:
            return self.app(environ, start_response)
        path = environ.get("PATH_INFO", "") or "/"
        if path.startswith(self.prefix):
            # Local/dev: request is /panviewer/search -> strip prefix, set SCRIPT_NAME
            environ["PATH_INFO"] = path[len(self.prefix) :] or "/"
            environ["SCRIPT_NAME"] = self.prefix
        else:
            # Production: proxy already stripped path; just set SCRIPT_NAME for url_for
            environ["SCRIPT_NAME"] = self.prefix
        return self.app(environ, start_response)

app.wsgi_app = PrefixMiddleware(app.wsgi_app, APPLICATION_ROOT)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_ROOT = os.path.join(BASE_DIR, "input")
DEFAULT_DATASET_ID = "wheat"

# Per-dataset config. Each dataset must have:
# - a SQLite DB created by `build_index.py`
# - optional per-cluster protein FASTA directory with files named like:
#   N0.HOG0000000.protein.fasta
DATASET_CONFIG = {
    "wheat": {
        # Pandagma pangenes mode (PANDAGMA pan IDs; no OrthoFinder-style OG/HOG hierarchy).
        # Sequences are fetched per gene_id from primary/<accession>.fa via read_fasta().
        "db_path": os.path.join(INPUT_ROOT, "wheat", "panwheat_pandagma.db"),
        "fasta_dir": None,
        "prot_dir": os.path.join(BASE_DIR, "primary", "prot"),
        "label": "Wheat",
        "mode": "pandagma",
        # Pandagma output does not have an OG→pangene hierarchy to drive the wheat OG picker UI.
        "enable_og_picker": False,
        # Disable extra OG/homeologue panels on the search results page (keep local synteny only).
        "enable_homeologue_panels": False,
    },
}


def load_external_db_links() -> list[dict[str, str]]:
    """
    Rows from static/links.csv for external DB buttons on the Genes tab.
    URLs may contain {gene_id} and/or {chr}, {start}, {end} placeholders.
    """
    path = os.path.join(BASE_DIR, "static", "links.csv")
    out: list[dict[str, str]] = []
    if not os.path.isfile(path):
        return out
    try:
        with open(path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                url = (row.get("url") or "").strip()
                if not url:
                    continue
                out.append(
                    {
                        "database": (row.get("database") or "").strip(),
                        "type": (row.get("type") or "").strip(),
                        "accession": (row.get("accession") or "").strip(),
                        "url": url,
                    }
                )
    except OSError:
        return out
    return out


EXTERNAL_DB_LINKS = load_external_db_links()


def load_about_stats(dataset_id: str) -> dict[str, str] | None:
    """
    Read input/<dataset>/dataset_stats.tsv (written by build_index or compute_dataset_stats).
    Returns display strings with thousands separators, or None if missing/invalid.
    """
    dataset_id = (dataset_id or "").strip().lower()
    path = os.path.join(INPUT_ROOT, dataset_id, "dataset_stats.tsv")
    if not os.path.isfile(path):
        return None
    try:
        with open(path, newline="") as f:
            row = next(csv.DictReader(f, delimiter="\t"), None)
        if not row:
            return None
        keys = ("accessions", "ogs", "hogs", "genes")
        if not all(k in row and str(row[k]).strip() for k in keys):
            return None
        return {k: f"{int(row[k]):,}" for k in keys}
    except (OSError, ValueError, TypeError):
        return None


def available_datasets():
    out = set()
    for ds_id, cfg in DATASET_CONFIG.items():
        if not os.path.exists(cfg.get("db_path", "")):
            continue

        mode = cfg.get("mode")
        if mode == "pandagma":
            # Sequences and synteny are served from the DB (e.g. protein_seqs, gene_coords).
            # primary/prot is only needed when building the index; do not require it at runtime.
            out.add(ds_id)
            continue

        # N0 + FASTA directory mode requires a per-cluster protein FASTA directory.
        fasta_dir = cfg.get("fasta_dir")
        if fasta_dir and os.path.isdir(fasta_dir):
            out.add(ds_id)
    return out


def get_dataset_config(dataset_id: str):
    dataset_id = (dataset_id or "").strip().lower()
    if dataset_id not in DATASET_CONFIG:
        dataset_id = DEFAULT_DATASET_ID
    return DATASET_CONFIG[dataset_id]


def get_db(dataset_id: str):
    cfg = get_dataset_config(dataset_id)
    conn = sqlite3.connect(cfg["db_path"])
    conn.row_factory = sqlite3.Row
    return conn


# gene_coords.chrom_index fast path for synteny (see build_index.assign_chrom_indices).
# sqlite3.Connection cannot be weak-referenced or given arbitrary attrs on some Python builds;
# cache by main database file path instead.
_gene_coords_meta_lock = threading.Lock()
_gene_coords_meta_by_db_path: dict[str, dict] = {}


def _main_sqlite_db_path(conn) -> str:
    """Absolute path of the main DB file, or ':memory:'."""
    rows = conn.execute("PRAGMA database_list").fetchall()
    for row in rows:
        if row[1] == "main":
            f = row[2]
            return os.path.abspath(f) if f else ":memory:"
    return ":memory:"


def _gene_coords_meta(conn) -> dict:
    """Cached PRAGMA + whether chrom_index is populated (for narrow synteny queries)."""
    key = _main_sqlite_db_path(conn)
    with _gene_coords_meta_lock:
        hit = _gene_coords_meta_by_db_path.get(key)
    if hit is not None:
        return hit

    cur = conn.cursor()
    cur.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='gene_coords'"
    )
    if not cur.fetchone():
        m = {
            "has_table": False,
            "window_ok": False,
            "select_list": "gene_id, accession, chr, start, end, strand, subgenome",
        }
        with _gene_coords_meta_lock:
            _gene_coords_meta_by_db_path[key] = m
        return m
    cur.execute("PRAGMA table_info(gene_coords)")
    cols = {r[1] for r in cur.fetchall()}
    has_ci = "chrom_index" in cols
    window_ok = False
    if has_ci:
        cur.execute(
            "SELECT 1 FROM gene_coords WHERE chrom_index IS NOT NULL LIMIT 1"
        )
        window_ok = cur.fetchone() is not None
    sel = "gene_id, accession, chr, start, end, strand, subgenome"
    if has_ci:
        sel += ", chrom_index"
    m = {"has_table": True, "window_ok": window_ok, "select_list": sel}
    with _gene_coords_meta_lock:
        _gene_coords_meta_by_db_path[key] = m
    return m


def extract_accession(gene_id):
    """Extract accession name from gene ID like HORVU.BONUS.PROJ.1HG00046750.1"""
    parts = gene_id.split(".")
    if len(parts) >= 3:
        return parts[1]
    return gene_id


def chinese_spring_gene_ids_from_rows(genes) -> list[str]:
    """All genes in this pangene cluster with Chinese Spring pangene accession (no space), case-insensitive, sorted."""
    seen: set[str] = set()
    out: list[str] = []
    for row in genes:
        acc = (row["accession"] or "").strip()
        if acc.lower() != "chinesespring":
            continue
        gid = (row["gene_id"] or "").strip()
        if gid and gid not in seen:
            seen.add(gid)
            out.append(gid)
    out.sort()
    return out


# Old name (singular); kept so stale bytecode / edits that still call it do not 500.
chinese_spring_gene_id_from_rows = chinese_spring_gene_ids_from_rows


def search_genes(query, dataset_id: str):
    """Search for genes matching the query (exact or partial)."""
    conn = get_db(dataset_id)
    cur = conn.cursor()
    cur.execute(
        "SELECT gene_id, hog, og, accession FROM genes WHERE gene_id = ? COLLATE NOCASE",
        (query.strip(),),
    )
    results = cur.fetchall()

    if not results:
        cur.execute(
            "SELECT gene_id, hog, og, accession FROM genes "
            "WHERE gene_id LIKE ? COLLATE NOCASE LIMIT 200",
            (f"%{query.strip()}%",),
        )
        results = cur.fetchall()

    conn.close()
    return results


def get_hog_genes(hog_id, dataset_id: str):
    """Get all genes belonging to a HOG (optional chr/start/end from gene_coords)."""
    conn = get_db(dataset_id)
    cur = conn.cursor()
    if table_exists(conn, "gene_coords"):
        cur.execute(
            """
            SELECT g.gene_id AS gene_id, g.accession AS accession,
                   gc.chr AS chr, gc.start AS start, gc.end AS end
            FROM genes g
            LEFT JOIN gene_coords gc
              ON g.gene_id = gc.gene_id AND g.accession = gc.accession
            WHERE g.hog = ?
            ORDER BY g.accession
            """,
            (hog_id,),
        )
    else:
        cur.execute(
            "SELECT gene_id, accession FROM genes WHERE hog = ? ORDER BY accession",
            (hog_id,),
        )
    results = cur.fetchall()
    conn.close()
    return results


def get_og_genes(og_id, dataset_id: str):
    """Get all genes belonging to an OG."""
    conn = get_db(dataset_id)
    cur = conn.cursor()
    cur.execute(
        "SELECT gene_id, hog, accession FROM genes WHERE og = ? ORDER BY hog, accession",
        (og_id,),
    )
    results = cur.fetchall()
    conn.close()
    return results


def build_wheat_og_picker_payload(
    conn,
    og_id: str,
    current_pangene_id: str,
    *,
    max_hogs: int = 200,
) -> dict | None:
    """
    Data for wheat pangene detail: all pangenes in the OG with genes + optional A/B/D dominants.
    Used for multi-select UI (accessions/genes/synteny scope).
    """
    cur = conn.cursor()
    cur.execute("SELECT COUNT(DISTINCT hog) AS n FROM genes WHERE og = ?", (og_id,))
    nh = int(cur.fetchone()["n"])
    if nh == 0:
        return None
    if nh > max_hogs:
        return {
            "og": og_id,
            "skipped": True,
            "reason": (
                f"This orthogroup has {nh} distinct pangenes; the picker supports up to {max_hogs}."
            ),
            "pangene_count": nh,
            "current_pangene": current_pangene_id,
        }

    coords_on = table_exists(conn, "gene_coords")
    cur.execute(
        "SELECT hog, gene_id, accession FROM genes WHERE og = ? ORDER BY hog, gene_id",
        (og_id,),
    )
    by_hog: dict[str, list[dict]] = defaultdict(list)
    for r in cur.fetchall():
        ge: dict = {
            "gene_id": r["gene_id"],
            "accession": r["accession"],
            "pangene": r["hog"],
            "og": og_id,
        }
        if coords_on:
            summ = _gene_coords_summary(cur, ge["gene_id"], ge["accession"])
            ge["coords"] = summ
            cr = fetch_coords_row(cur, ge["gene_id"], ge["accession"])
            if cr:
                sg_raw = (cr.get("subgenome") or "").strip().upper()
                ge["subgenome"] = sg_raw or None
                ge["chr"] = cr["chr"]
                ge["start"] = int(cr["start"])
                ge["end"] = int(cr["end"])
                st = (cr.get("strand") or "").strip()
                ge["strand"] = st if st else None
            else:
                ge["subgenome"] = None
                ge["chr"] = ge["start"] = ge["end"] = None
                ge["strand"] = None
        else:
            ge["coords"] = None
            ge["subgenome"] = None
            ge["chr"] = ge["start"] = ge["end"] = ge["strand"] = None
        by_hog[r["hog"]].append(ge)
    hogs_sorted = sorted(by_hog.keys())

    panel = None
    if coords_on:
        panel = build_wheat_og_pangene_triad_panel(
            conn, og_id, current_pangene_id, max_hogs=max_hogs
        )
    row_by_pangene: dict[str, dict] = {}
    if panel and not panel.get("skipped"):
        for r in panel["rows"]:
            row_by_pangene[r["pangene"]] = r

    items: list[dict] = []
    for h in hogs_sorted:
        pr = row_by_pangene.get(h)
        gl = by_hog[h]
        items.append(
            {
                "pangene": h,
                "is_current": h == current_pangene_id,
                "dominant": pr["dominant"] if pr else "?",
                "gene_count": pr["gene_count"] if pr else len(gl),
                "counts": pr["counts"]
                if pr
                else {"A": 0, "B": 0, "D": 0, "?": len(gl)},
                "genes": gl,
            }
        )

    items.sort(
        key=lambda x: (0 if x["pangene"] == current_pangene_id else 1, x["pangene"])
    )

    return {
        "og": og_id,
        "skipped": False,
        "pangene_items": items,
        "current_pangene": current_pangene_id,
    }


def get_hog_info(hog_id, dataset_id: str):
    """Get pangene cluster metadata (SQLite table hog_info; column hog = pan / cluster id)."""
    conn = get_db(dataset_id)
    cur = conn.cursor()
    cur.execute("SELECT * FROM hog_info WHERE hog = ?", (hog_id,))
    result = cur.fetchone()
    conn.close()
    return result


_PROT_INDEX_CACHE: dict[str, Any] = {}
_PROT_INDEX_CACHE_LOCK = threading.Lock()


def _get_protein_index(prot_fasta_path: str):
    """Cache Biopython SeqIO.index for random gene_id access."""
    with _PROT_INDEX_CACHE_LOCK:
        hit = _PROT_INDEX_CACHE.get(prot_fasta_path)
        if hit is not None:
            return hit
    # SeqIO.index builds an offset map in memory; sequences are loaded lazily.
    # In restricted environments, SeqIO.index can fail (e.g., OpenMP/shared-memory).
    # Fallback to streaming parse into a dict so pandagma mode still works.
    try:
        idx = SeqIO.index(prot_fasta_path, "fasta")
    except Exception:
        seqs: dict[str, str] = {}
        for record in SeqIO.parse(prot_fasta_path, "fasta"):
            seqs[record.id] = str(record.seq)
        idx = seqs
    with _PROT_INDEX_CACHE_LOCK:
        _PROT_INDEX_CACHE[prot_fasta_path] = idx
    return idx


def read_fasta(hog_id, dataset_id: str):
    """
    Read protein sequences for a given group id (Pandagma pan_id in Pandagma mode).
    """
    cfg = get_dataset_config(dataset_id)

    fasta_dir = cfg.get("fasta_dir")
    fasta_path = (
        os.path.join(fasta_dir, f"{hog_id}.protein.fasta")
        if fasta_dir
        else None
    )
    if fasta_path and os.path.exists(fasta_path):
        sequences = {}
        for record in SeqIO.parse(fasta_path, "fasta"):
            sequences[record.id] = str(record.seq)
        return sequences

    if cfg.get("mode") == "pandagma":
        conn = get_db(dataset_id)
        try:
            if not table_exists(conn, "protein_seqs"):
                return {}
            cur = conn.cursor()
            cur.execute(
                """
                SELECT g.gene_id, p.seq_comp
                FROM genes g
                JOIN protein_seqs p ON p.gene_id = g.gene_id
                WHERE g.hog = ?
                ORDER BY g.accession, g.gene_id
                """,
                (hog_id,),
            )
            rows = cur.fetchall()
        finally:
            conn.close()

        sequences: dict[str, str] = {}
        for r in rows:
            gid = (r["gene_id"] or "").strip()
            comp = r["seq_comp"]
            if not gid or comp is None:
                continue
            try:
                sequences[gid] = zlib.decompress(comp).decode("utf-8")
            except Exception:
                continue

        return sequences

    return {}


def collect_sequences_for_gene_ids(
    conn: sqlite3.Connection, gene_ids: list[str], dataset_id: str
) -> tuple[dict[str, str], dict[str, str]]:
    """
    Load ungapped protein sequences for the given gene IDs (any hog/pan in DB).
    Returns (gene_id -> sequence, gene_id -> accession).
    """
    cur = conn.cursor()
    hog_to_ids: dict[str, set[str]] = {}
    acc_map: dict[str, str] = {}
    for gid in gene_ids:
        row = None
        for vid in gene_id_coord_variants(gid):
            cur.execute(
                """
                SELECT gene_id, hog, accession FROM genes
                WHERE gene_id = ? COLLATE NOCASE LIMIT 1
                """,
                (vid,),
            )
            r = cur.fetchone()
            if r:
                row = r
                break
        if not row:
            continue
        db_g = (row["gene_id"] or "").strip()
        hog = (row["hog"] or "").strip()
        acc = (row["accession"] or "").strip()
        if not db_g or not hog:
            continue
        hog_to_ids.setdefault(hog, set()).add(db_g)
        acc_map[db_g] = acc or db_g

    merged: dict[str, str] = {}
    for hog, idset in hog_to_ids.items():
        seqs = read_fasta(hog, dataset_id)
        if not seqs:
            continue
        low = {k.lower(): k for k in seqs}
        for db_g in idset:
            if db_g in seqs:
                merged[db_g] = str(seqs[db_g]).strip()
                continue
            kk = low.get(db_g.lower())
            if kk:
                merged[db_g] = str(seqs[kk]).strip()
    return merged, acc_map


# FAMSA2: -t 0 = half of logical cores (see FAMSA README).
FAMSA_TIMEOUT_SEC = 600


def run_famsa(sequences, *, timeout_sec=FAMSA_TIMEOUT_SEC):
    """Align with FAMSA only using an NJ guide tree.

    Exports NJ Newick (-gt nj -gt_export), then aligns with -gt import so the
    tree matches the progressive alignment. Returns (aligned_dict, newick_str).
    On missing binary, failure, or timeout, returns (original sequences, "").
    """
    log = logging.getLogger(__name__)
    if len(sequences) <= 1:
        return sequences, ""

    famsa_bin = shutil.which("famsa")
    if not famsa_bin:
        log.warning("famsa not found in PATH")
        return sequences, ""

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    ) as tmp_in:
        for sid, seq in sequences.items():
            clean_seq = str(seq).rstrip("*")
            tmp_in.write(f">{sid}\n{clean_seq}\n")
        tmp_in_path = tmp_in.name

    tmp_dnd = tmp_in_path + ".nj.dnd"
    tmp_aln = tmp_in_path + ".famsa.aln"

    def _cleanup():
        for p in (tmp_in_path, tmp_dnd, tmp_aln):
            if os.path.exists(p):
                try:
                    os.unlink(p)
                except OSError:
                    pass

    # Do not pass -keep_duplicates: FAMSA 1.x ignores unknown flags and treats the
    # next token as the input path, yielding "Unable to open input file -keep_duplicates".
    famsa_base = [famsa_bin, "-t", "0"]
    export_cmd = famsa_base + ["-gt", "nj", "-gt_export", tmp_in_path, tmp_dnd]
    import_cmd = famsa_base + ["-gt", "import", tmp_dnd, tmp_in_path, tmp_aln]

    try:
        try:
            r1 = subprocess.run(
                export_cmd,
                capture_output=True,
                text=True,
                timeout=timeout_sec,
            )
        except FileNotFoundError:
            log.warning("famsa binary disappeared from PATH")
            return sequences, ""

        if r1.returncode != 0 or not os.path.isfile(tmp_dnd):
            err = (r1.stderr or r1.stdout or "").strip()[:800]
            log.warning("famsa NJ export failed: %s", err or "(no output)")
            return sequences, ""

        with open(tmp_dnd, encoding="utf-8", errors="replace") as f:
            newick = f.read().strip()

        try:
            r2 = subprocess.run(
                import_cmd,
                capture_output=True,
                text=True,
                timeout=timeout_sec,
            )
        except FileNotFoundError:
            return sequences, ""

        if r2.returncode != 0 or not os.path.isfile(tmp_aln):
            err = (r2.stderr or r2.stdout or "").strip()[:800]
            log.warning("famsa import alignment failed: %s", err or "(no output)")
            return sequences, ""

        aligned = {}
        for record in SeqIO.parse(tmp_aln, "fasta"):
            aligned[record.id] = str(record.seq)
        if not aligned:
            return sequences, ""
        return aligned, newick
    except subprocess.TimeoutExpired:
        log.warning("famsa timed out after %s s", timeout_sec)
        return sequences, ""
    finally:
        _cleanup()


def compute_consensus(aligned_sequences):
    """Compute consensus sequence from aligned sequences."""
    if not aligned_sequences:
        return ""
    seqs = list(aligned_sequences.values())
    aln_len = max(len(s) for s in seqs)
    consensus = []
    for i in range(aln_len):
        freq = defaultdict(int)
        for s in seqs:
            aa = s[i] if i < len(s) else "-"
            if aa != "-":
                freq[aa] += 1
        consensus.append(max(freq, key=freq.get) if freq else "-")
    return "".join(consensus)


def get_tree_leaf_order(newick_str):
    """Parse a Newick tree and return leaf names in display order (top to bottom)."""
    if not newick_str:
        return []
    try:
        tree = Phylo.read(io.StringIO(newick_str), "newick")
        return [clade.name for clade in tree.get_terminals()]
    except Exception:
        return []


def compute_conservation(aligned_sequences):
    """
    Per-column conservation from alignment (no R/bio3d).
    Score = fraction of most common non-gap residue at each column (0-1, 1 = fully conserved).
    """
    if not aligned_sequences:
        return []
    seqs = list(aligned_sequences.values())
    aln_len = max(len(s) for s in seqs)
    out = []
    for i in range(aln_len):
        freq = defaultdict(int)
        for s in seqs:
            aa = s[i] if i < len(s) else "-"
            if aa != "-":
                freq[aa] += 1
        n = sum(freq.values())
        if n == 0:
            out.append(0.0)
        else:
            out.append(max(freq.values()) / n)
    return out


def _clean_dna(s: str) -> str:
    return "".join(b for b in (s or "").upper() if b in "ATCG")


def _translate_codon(codon: str) -> str | None:
    c = codon.upper()
    if len(c) != 3 or set(c) - set("ATCG"):
        return None
    try:
        aa = str(Seq(c).translate(table=1))
    except Exception:
        return None
    if aa == "*":
        return None
    return aa


def _ng86_sites_per_codon(codon: str) -> tuple[float, float]:
    """Nei–Gojobori-style expected synonymous / nonsynonymous sites for one codon."""
    c = codon.upper()
    if len(c) != 3:
        return 0.0, 0.0
    aa = _translate_codon(c)
    if aa is None:
        return 0.0, 0.0
    s_frac = 0.0
    for i in range(3):
        local_syn = 0
        for b in "ATGC":
            if b == c[i]:
                continue
            mut = c[:i] + b + c[i+1:]
            aa2 = _translate_codon(mut)
            if aa2 is not None and aa2 == aa:
                local_syn += 1
        s_frac += local_syn / 3.0
    n_frac = 3.0 - s_frac
    return s_frac, n_frac


def map_protein_alignment_to_codons(aln_prot: str, cds: str) -> list[str] | None:
    """One codon string per alignment column (including '---' for AA gaps)."""
    cds_clean = _clean_dna(cds)
    n_aa = sum(1 for x in aln_prot if x != "-")
    if n_aa == 0:
        return None
    need = 3 * n_aa
    if len(cds_clean) < need:
        return None
    cds_use = cds_clean[:need]
    out: list[str] = []
    ci = 0
    for aa in aln_prot:
        if aa == "-":
            out.append("---")
        else:
            out.append(cds_use[ci : ci + 3])
            ci += 3
    return out


def _ng86_column_contrib(cr: str, cq: str) -> tuple[float, float, int, int]:
    """Per-column (S_sites, N_sites, Sd, Nd) for Nei–Gojobori accumulation."""
    if cr == "---" or cq == "---":
        return 0.0, 0.0, 0, 0
    if len(cr) != 3 or len(cq) != 3:
        return 0.0, 0.0, 0, 0
    if "N" in cr.upper() or "N" in cq.upper():
        return 0.0, 0.0, 0, 0
    if _translate_codon(cr) is None or _translate_codon(cq) is None:
        return 0.0, 0.0, 0, 0
    sr, nr = _ng86_sites_per_codon(cr)
    sq, nq = _ng86_sites_per_codon(cq)
    S_sites = (sr + sq) / 2.0
    N_sites = (nr + nq) / 2.0
    if cr == cq:
        return S_sites, N_sites, 0, 0
    tr, tq = _translate_codon(cr), _translate_codon(cq)
    if tr is None or tq is None:
        return 0.0, 0.0, 0, 0
    if tr == tq:
        return S_sites, N_sites, 1, 0
    return S_sites, N_sites, 0, 1


def _omega_from_ng86_totals(
    Nd: int, Sd: int, N_sites: float, S_sites: float
) -> float | None:
    if S_sites <= 0 and N_sites <= 0:
        return None
    Ka = Nd / N_sites if N_sites > 0 else None
    Ks = Sd / S_sites if S_sites > 0 else None
    if Ka is not None and Ks is not None and Ks > 0:
        return float(Ka / Ks)
    return None


def _ng86_prefix_sums(
    ref_codons: list[str], q_codons: list[str]
) -> tuple[list[float], list[float], list[int], list[int]]:
    """Prefix sums so any window [lo, hi] inclusive is O(1)."""
    n = len(ref_codons)
    ps_s = [0.0] * (n + 1)
    ps_n = [0.0] * (n + 1)
    ps_sd = [0] * (n + 1)
    ps_nd = [0] * (n + 1)
    for i in range(n):
        ss, ns, sd, nd = _ng86_column_contrib(ref_codons[i], q_codons[i])
        ps_s[i + 1] = ps_s[i] + ss
        ps_n[i + 1] = ps_n[i] + ns
        ps_sd[i + 1] = ps_sd[i] + sd
        ps_nd[i + 1] = ps_nd[i] + nd
    return ps_s, ps_n, ps_sd, ps_nd


def pairwise_kaks_ng86(ref_codons: list[str], q_codons: list[str]) -> dict[str, Any]:
    """Pairwise Nei–Gojobori–style Ka, Ks, omega from aligned codon rows."""
    if len(ref_codons) != len(q_codons):
        return {"ok": False, "error": "length_mismatch"}
    ps_s, ps_n, ps_sd, ps_nd = _ng86_prefix_sums(ref_codons, q_codons)
    n = len(ref_codons)
    S_sites = ps_s[n]
    N_sites = ps_n[n]
    Sd = ps_sd[n]
    Nd = ps_nd[n]
    if S_sites <= 0 and N_sites <= 0:
        return {"ok": False, "error": "no_valid_sites"}
    Ka = Nd / N_sites if N_sites > 0 else None
    Ks = Sd / S_sites if S_sites > 0 else None
    omega: float | None = None
    if Ka is not None and Ks is not None and Ks > 0:
        omega = Ka / Ks
    return {
        "ok": True,
        "Nd": Nd,
        "Sd": Sd,
        "N_sites": N_sites,
        "S_sites": S_sites,
        "Ka": Ka,
        "Ks": Ks,
        "omega": omega,
    }


def load_cds_for_gene(cur, gene_id: str) -> str | None:
    """Decompress CDS from cds_seqs, trying transcript ID variants."""
    for vid in gene_id_coord_variants(gene_id):
        cur.execute(
            "SELECT seq_comp FROM cds_seqs WHERE gene_id = ? LIMIT 1",
            (vid,),
        )
        r = cur.fetchone()
        if r and r["seq_comp"]:
            try:
                return zlib.decompress(r["seq_comp"]).decode("utf-8")
            except Exception:
                continue
    return None


def compute_kaks_vs_reference(
    conn,
    hog_id: str,
    ref_gene_id: str,
    aligned_proteins: dict[str, str],
) -> list[dict[str, Any]]:
    """
    Ka/Ks for each non-reference gene vs ref using CDS mapped onto the given protein alignment.
    aligned_proteins must match the in-app MSA (same gaps as FAMSA output).
    """
    if not table_exists(conn, "cds_seqs"):
        return []
    cur = conn.cursor()
    cur.execute("SELECT gene_id FROM genes WHERE hog = ?", (hog_id,))
    allowed = {r["gene_id"] for r in cur.fetchall()}
    ref_aln = aligned_proteins.get(ref_gene_id)
    if not ref_aln:
        return []
    ref_cds = load_cds_for_gene(cur, ref_gene_id)
    if not ref_cds:
        return []
    ref_codons = map_protein_alignment_to_codons(ref_aln, ref_cds)
    if not ref_codons:
        return []

    rows: list[dict[str, Any]] = []
    for gid, aln in sorted(aligned_proteins.items()):
        if gid == ref_gene_id:
            continue
        if gid not in allowed:
            continue
        cds = load_cds_for_gene(cur, gid)
        if not cds:
            rows.append(
                {
                    "gene_id": gid,
                    "ok": False,
                    "error": "no_cds",
                }
            )
            continue
        q_codons = map_protein_alignment_to_codons(aln, cds)
        if not q_codons:
            rows.append(
                {
                    "gene_id": gid,
                    "ok": False,
                    "error": "cds_length_mismatch",
                }
            )
            continue
        stats = pairwise_kaks_ng86(ref_codons, q_codons)
        if not stats.get("ok"):
            rows.append(
                {
                    "gene_id": gid,
                    "ok": False,
                    "error": stats.get("error", "compute_failed"),
                }
            )
            continue
        rows.append(
            {
                "gene_id": gid,
                "ok": True,
                "Ka": stats["Ka"],
                "Ks": stats["Ks"],
                "omega": stats["omega"],
                "Nd": stats["Nd"],
                "Sd": stats["Sd"],
                "N_sites": stats["N_sites"],
                "S_sites": stats["S_sites"],
            }
        )
    return rows


def compute_sliding_kaks_profile(
    conn,
    hog_id: str,
    ref_gene_id: str,
    aligned_proteins: dict[str, str],
    *,
    window_half: int = 5,
) -> list[float | None]:
    """
    Per alignment column: mean pairwise ω (Nei–Gojobori) in a sliding window of codons
    (ref vs each other active gene with CDS). Approximates where constraint vs
    diversification differs along the protein; not a formal site model.
    """
    if not table_exists(conn, "cds_seqs"):
        return []
    cur = conn.cursor()
    cur.execute("SELECT gene_id FROM genes WHERE hog = ?", (hog_id,))
    allowed = {r["gene_id"] for r in cur.fetchall()}
    ref_aln = aligned_proteins.get(ref_gene_id)
    if not ref_aln:
        return []
    alen = len(ref_aln)
    if ref_gene_id not in allowed:
        return [None] * alen

    ref_cds = load_cds_for_gene(cur, ref_gene_id)
    if not ref_cds:
        return [None] * alen
    ref_codons = map_protein_alignment_to_codons(ref_aln, ref_cds)
    if not ref_codons or len(ref_codons) != alen:
        return [None] * alen

    others = [
        g
        for g in aligned_proteins
        if g != ref_gene_id and g in allowed and isinstance(aligned_proteins.get(g), str)
    ]
    q_codons_map: dict[str, list[str]] = {}
    for gid in others:
        cds = load_cds_for_gene(cur, gid)
        if not cds:
            continue
        qc = map_protein_alignment_to_codons(aligned_proteins[gid], cds)
        if qc and len(qc) == alen:
            q_codons_map[gid] = qc

    if not q_codons_map:
        return [None] * alen

    wh = max(1, min(int(window_half), 80))

    # Prefix sums per (ref vs query): O(alen) each; sliding windows O(1) — not O(window) per column.
    pair_prefix: dict[str, tuple[list[float], list[float], list[int], list[int]]] = {}
    for gid, qc in q_codons_map.items():
        pair_prefix[gid] = _ng86_prefix_sums(ref_codons, qc)

    profile: list[float | None] = []
    for i in range(alen):
        if ref_codons[i] == "---":
            profile.append(None)
            continue
        lo = max(0, i - wh)
        hi = min(alen - 1, i + wh)
        omegas: list[float] = []
        for _gid, pfx in pair_prefix.items():
            ps_s, ps_n, ps_sd, ps_nd = pfx
            S_sites = ps_s[hi + 1] - ps_s[lo]
            N_sites = ps_n[hi + 1] - ps_n[lo]
            Sd = ps_sd[hi + 1] - ps_sd[lo]
            Nd = ps_nd[hi + 1] - ps_nd[lo]
            om = _omega_from_ng86_totals(Nd, Sd, N_sites, S_sites)
            if om is not None and math.isfinite(om):
                omegas.append(om)
        if not omegas:
            profile.append(None)
        else:
            profile.append(sum(omegas) / len(omegas))
    return profile


# --- Wheat coords / homeologues / synteny (optional `gene_coords` table) ---


def table_exists(conn, name: str) -> bool:
    cur = conn.cursor()
    cur.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=? LIMIT 1",
        (name,),
    )
    return cur.fetchone() is not None


def _ungapped_protein_len(seq: str) -> int:
    return len(re.sub(r"[.\-]", "", str(seq or "")))


def _pad_or_trim_porter_ss(ss: str | None, n: int) -> str:
    """Match Porter string length to ungapped protein length (pad with '-', trim excess)."""
    s = (ss or "").replace(" ", "")
    if n <= 0:
        return ""
    if len(s) >= n:
        return s[:n]
    return s + "-" * (n - len(s))


def fetch_porter6_raw_for_genes(
    conn: sqlite3.Connection, gene_ids: list[str], seq_lens: dict[str, int]
) -> dict[str, dict[str, str]]:
    """
    Per-gene Porter6 q3/q8 strings (same length as ungapped protein).
    Missing DB rows or missing table → all '-' for that length.

    Uses batched ``WHERE gene_id IN (...)`` so SQLite can use ``idx_porter6_ss_gene_id``.
    (Per-row ``WHERE gene_id = ? COLLATE NOCASE`` cannot use that index and devolves to
    full table scans on multi-million-row ``porter6_ss`` tables.)
    """
    out: dict[str, dict[str, str]] = {}
    if not gene_ids:
        return out
    has_table = table_exists(conn, "porter6_ss")
    cur = conn.cursor()

    # Unique IDs preserving order (alignment keys)
    seen: set[str] = set()
    unique_ids: list[str] = []
    for gid in gene_ids:
        if gid not in seen:
            seen.add(gid)
            unique_ids.append(gid)

    qmap: dict[str, tuple[str, str]] = {}
    qmap_lower: dict[str, tuple[str, str]] = {}
    if has_table and unique_ids:
        chunk_size = 400
        for i in range(0, len(unique_ids), chunk_size):
            chunk = unique_ids[i : i + chunk_size]
            ph = ",".join("?" * len(chunk))
            cur.execute(
                f"SELECT gene_id, q3, q8 FROM porter6_ss WHERE gene_id IN ({ph})",
                chunk,
            )
            for r in cur.fetchall():
                g = (r["gene_id"] or "").strip()
                if not g:
                    continue
                tup = ((r["q3"] or ""), (r["q8"] or ""))
                qmap[g] = tup
                qmap_lower.setdefault(g.lower(), tup)

    for gid in gene_ids:
        n = max(1, int(seq_lens.get(gid, 1)))
        if not has_table:
            out[gid] = {"q3": "-" * n, "q8": "-" * n}
            continue
        row = qmap.get(gid) or qmap_lower.get(gid.lower())
        if row:
            q3 = _pad_or_trim_porter_ss(row[0], n)
            q8 = _pad_or_trim_porter_ss(row[1], n)
        else:
            q3 = q8 = "-" * n
        out[gid] = {"q3": q3, "q8": q8}
    return out


def porter6_raw_for_aligned_sequences(
    conn: sqlite3.Connection | None, aligned: dict[str, str]
) -> dict[str, dict[str, str]]:
    """Build porter6 payload for keys present in an alignment dict."""
    if not aligned:
        return {}
    seq_lens = {g: _ungapped_protein_len(s) for g, s in aligned.items()}
    if conn is None:
        return {
            g: {"q3": "-" * max(1, seq_lens.get(g, 1)), "q8": "-" * max(1, seq_lens.get(g, 1))}
            for g in aligned
        }
    return fetch_porter6_raw_for_genes(conn, list(aligned.keys()), seq_lens)


def porter_q3_to_secstructartist_hsl(s: str) -> str:
    """Map Porter 3-state (H/E/C/…) to secstructartist HSL alphabet (H/S/L)."""
    out: list[str] = []
    for c in s:
        u = c.upper()
        if u == "H":
            out.append("H")
        elif u == "E":
            out.append("S")
        elif u in ("C", "-", ".", " "):
            out.append("L")
        else:
            out.append(u)
    return "".join(out)


# Element colors for secstructartist (Conservation /api/ss_plot). Kept in sync with
# templates/pangene_detail.html ssColorForChar (alignment + conservation strip).
_SS_ARTIST_HEIGHT = 0.82
_SS_EDGE = "#555555"


def secstructartist_svg_from_rows(
    rows: list[str], mode: str
) -> tuple[str | None, str | None]:
    """
    Render a multi-row secondary-structure diagram as SVG (in-memory; no temp files).
    Q3 uses SecStructArtist (HSL); Q8 uses SecStructArtistDssp (DSSP letters).
    Returns (svg_xml, error_message).
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from secstructartist.interface.artist import SecStructArtist
        from secstructartist.interface.artist_dssp import SecStructArtistDssp
    except ImportError as e:
        return None, f"import_error:{e!s}"

    if not rows:
        return None, "no_rows"
    if mode not in ("q3", "q8"):
        return None, "bad_mode"

    # Pad every row to the same length (50 residues / row from the client) so the last
    # row is not stretched across the full axis; trailing positions are blank.
    pad_w = max(len(r) for r in rows)
    rows = [
        (r if len(r) >= pad_w else r + (" " * (pad_w - len(r))))[:pad_w] for r in rows
    ]

    n = len(rows)
    maxlen = pad_w
    fig_w = max(5.0, 0.11 * maxlen)
    # Slightly shorter rows than default (was 0.52 * n, min 0.9) for a less tall diagram.
    fig_h = max(0.72, 0.42 * n)
    fig, axes = plt.subplots(n, 1, figsize=(fig_w, fig_h))
    if n == 1:
        axes = [axes]

    if mode == "q3":
        artist = SecStructArtist(
            height=_SS_ARTIST_HEIGHT,
            helix_kwargs={
                "fillcolor": "#fae8e6",
                "shadecolor": "#f5b7b1",
                "linecolor": _SS_EDGE,
            },
            sheet_kwargs={"fillcolor": "#aed6f1", "linecolor": _SS_EDGE},
            loop_kwargs={"linecolor": "#868e96"},
        )
        for ax, row in zip(axes, rows):
            mapped = porter_q3_to_secstructartist_hsl(row)
            artist.draw(mapped, None, 0.0, ax=ax)
    else:
        artist = SecStructArtistDssp(
            height=_SS_ARTIST_HEIGHT,
            H_kwargs={
                "fillcolor": "#fae8e6",
                "shadecolor": "#f5b7b1",
                "linecolor": _SS_EDGE,
            },
            G_kwargs={
                "fillcolor": "#fef5e4",
                "shadecolor": "#fad7a0",
                "linecolor": _SS_EDGE,
            },
            I_kwargs={
                "fillcolor": "#f0e8f6",
                "shadecolor": "#d7bde2",
                "linecolor": _SS_EDGE,
            },
            E_kwargs={"fillcolor": "#aed6f1", "linecolor": _SS_EDGE},
            B_kwargs={"fillcolor": "#a3e4d7", "linecolor": _SS_EDGE},
            S_kwargs={"linecolor": "#5d6d7e"},
            T_kwargs={"linecolor": "#2980b9"},
            C_Kwargs={"linecolor": "#868e96"},
        )
        for ax, row in zip(axes, rows):
            mapped = "".join(
                "-" if c.isspace() else (c.upper() if c not in "-." else "-")
                for c in row
            )
            artist.draw(mapped, None, 0.0, ax=ax)

    # Y: no ticks/labels; X: ticks/labels on bottom subplot only.
    for i, ax in enumerate(axes):
        ax.set_ylabel("")
        ax.set_yticks([])
        ax.tick_params(axis="y", left=False, labelleft=False, right=False, labelright=False)
        if i < n - 1:
            ax.tick_params(axis="x", bottom=False, labelbottom=False)
            ax.set_xticklabels([])

    plt.subplots_adjust(hspace=0.36, left=0.03, right=0.97, top=0.96, bottom=0.04)
    buf = io.BytesIO()
    fig.savefig(buf, format="svg", bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)
    return buf.getvalue().decode("utf-8"), None


def gene_id_coord_variants(gene_id: str) -> list[str]:
    """Match coords TSV gene_id (often no transcript suffix) to genes table IDs."""
    g = (gene_id or "").strip()
    if not g:
        return []
    variants = [g]
    if re.search(r"\.\d+$", g):
        variants.append(re.sub(r"\.\d+$", "", g))
    else:
        variants.append(f"{g}.1")
    out, seen = [], set()
    for v in variants:
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


def _norm_gid(g: str) -> str:
    return re.sub(r"\.\d+$", "", (g or "").strip()).lower()


def build_gene_accession_lookup(hog_id: str, dataset_id: str) -> dict[str, str]:
    """
    Map alignment/FASTA IDs (incl. transcript variants) to genes.accession.
    Wheat IDs are not HORVU-style dotted names; client must not infer accession from split('.')[1].
    """
    m: dict[str, str] = {}
    for g in get_hog_genes(hog_id, dataset_id):
        acc = g["accession"]
        gid = g["gene_id"]
        for v in gene_id_coord_variants(gid):
            m[v] = acc
    return m


def build_multi_hog_gene_accession_lookup(
    hog_ids: list[str], dataset_id: str
) -> dict[str, str]:
    """Union of gene_id → accession for many HOGs (for merged alignment view)."""
    m: dict[str, str] = {}
    conn = get_db(dataset_id)
    try:
        cur = conn.cursor()
        for hid in hog_ids:
            cur.execute(
                "SELECT gene_id, accession FROM genes WHERE hog = ?",
                (hid,),
            )
            for r in cur.fetchall():
                for v in gene_id_coord_variants(r["gene_id"]):
                    m[v] = r["accession"]
    finally:
        conn.close()
    return m


def fetch_gene_ortho_batch(cur, gene_ids: list[str]) -> dict[str, dict | None]:
    """Orthology row per raw gene_id; one SQL query for the whole neighborhood."""
    if not gene_ids:
        return {}
    variants: list[str] = []
    seen_v: set[str] = set()
    for gid in gene_ids:
        for v in gene_id_coord_variants(gid):
            if v not in seen_v:
                seen_v.add(v)
                variants.append(v)
    if not variants:
        return {g: None for g in gene_ids}
    ph = ",".join("?" * len(variants))
    cur.execute(
        f"SELECT gene_id, hog, og, accession FROM genes WHERE gene_id IN ({ph})",
        variants,
    )
    by_norm: dict[str, dict] = {}
    for r in cur.fetchall():
        nk = _norm_gid(r["gene_id"])
        if nk not in by_norm:
            by_norm[nk] = dict(r)
    out: dict[str, dict | None] = {}
    for gid in gene_ids:
        ortho = None
        for v in gene_id_coord_variants(gid):
            row = by_norm.get(_norm_gid(v))
            if row:
                ortho = row
                break
        out[gid] = ortho
    return out


def fetch_coords_row(cur, gene_id: str, accession: str):
    meta = _gene_coords_meta(cur.connection)
    if not meta["has_table"]:
        return None
    sel = meta["select_list"]
    for vid in gene_id_coord_variants(gene_id):
        cur.execute(
            f"""
            SELECT {sel}
            FROM gene_coords
            WHERE gene_id = ? AND accession = ? COLLATE NOCASE
            """,
            (vid, accession),
        )
        r = cur.fetchone()
        if r:
            return dict(r)
    return None


def _norm_acc_key(acc: str) -> str:
    return (acc or "").strip().lower()


def fetch_coords_rows_batch(
    cur, gene_acc_pairs: list[tuple[str, str]]
) -> dict[str, dict | None]:
    """
    Like fetch_coords_row for many (gene_id, accession) pairs: one or few queries
    instead of O(pairs) round-trips. Keys in the result are the input gene_id strings.
    """
    meta = _gene_coords_meta(cur.connection)
    if not meta["has_table"]:
        return {g: None for g, _ in gene_acc_pairs}
    sel = meta["select_list"]

    uniq_pairs: list[tuple[str, str]] = []
    seen_q: set[tuple[str, str]] = set()
    for gid, acc in gene_acc_pairs:
        for vid in gene_id_coord_variants(gid):
            t = (vid, acc)
            if t not in seen_q:
                seen_q.add(t)
                uniq_pairs.append(t)

    fetched: dict[tuple[str, str], dict] = {}
    chunk_size = 60
    for i in range(0, len(uniq_pairs), chunk_size):
        chunk = uniq_pairs[i : i + chunk_size]
        parts: list[str] = []
        params: list[str] = []
        for vid, acc in chunk:
            parts.append("(gene_id = ? AND accession = ? COLLATE NOCASE)")
            params.extend([vid, acc])
        cur.execute(
            f"SELECT {sel} FROM gene_coords WHERE {' OR '.join(parts)}",
            params,
        )
        for r in cur.fetchall():
            rd = dict(r)
            gk = (rd["gene_id"], _norm_acc_key(rd["accession"]))
            fetched[gk] = rd

    out: dict[str, dict | None] = {}
    for gid, acc in gene_acc_pairs:
        cr = None
        for vid in gene_id_coord_variants(gid):
            cr = fetched.get((vid, _norm_acc_key(acc)))
            if cr:
                break
        out[gid] = cr
    return out


def _gene_coords_summary(cur, gene_id: str, accession: str) -> str | None:
    """Single-line chr:start-end [strand] for UI, or None."""
    c = fetch_coords_row(cur, gene_id, accession)
    if not c:
        return None
    frag = f"{c['chr']}:{int(c['start'])}-{int(c['end'])}"
    st = (c.get("strand") or "").strip()
    if st:
        frag += f" {st}"
    return frag


def fetch_gene_ortho(cur, gene_id: str):
    for vid in gene_id_coord_variants(gene_id):
        cur.execute(
            """
            SELECT gene_id, hog, og, accession
            FROM genes WHERE gene_id = ? COLLATE NOCASE
            """,
            (vid,),
        )
        r = cur.fetchone()
        if r:
            return dict(r)
    return None


def _synteny_color_for_key(key: str) -> str:
    """Stable HSL fill for a pangene (pan) or OG id.

    Colors do not depend on which other groups appear in the same synteny window,
    so the same pan / type looks identical across stacked rows (e.g. per-gene tracks).
    """
    k = (key or "").strip()
    if not k or k == "—":
        return "#c8c8c8"
    # Discrete palette: more distinct swatches than raw HSL jitter.
    # Chosen to keep decent contrast against the neutral grays used for missing/filtered rows.
    palette = [
        "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666",
        "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
        "#bcbd22", "#17becf", "#393b79", "#637939", "#843c39", "#8c6d31", "#ad494a", "#f7b6d2",
        "#c7e9c0", "#de9f76", "#bcbddc", "#9ecae1", "#fdd0a2", "#f768a1",
    ]
    digest = hashlib.md5(k.encode("utf-8"), usedforsecurity=False).digest()
    idx = ((digest[0] << 8) | digest[1]) % len(palette)
    return palette[idx]


def build_wheat_og_subgenome_table(conn, og_id: str, *, limit: int = 600) -> dict:
    """Partition genes in the same OG by subgenome (coords); may truncate."""
    cur = conn.cursor()
    cur.execute("SELECT COUNT(*) AS c FROM genes WHERE og = ?", (og_id,))
    total = int(cur.fetchone()["c"])
    cur.execute(
        """
        SELECT gene_id, hog, og, accession
        FROM genes WHERE og = ?
        ORDER BY accession, gene_id
        LIMIT ?
        """,
        (og_id, limit),
    )
    og_genes = [dict(r) for r in cur.fetchall()]
    by_sg: dict[str, list[dict]] = {"A": [], "B": [], "D": []}
    missing_coords: list[str] = []
    unplaced: list[dict] = []

    for g in og_genes:
        c = fetch_coords_row(cur, g["gene_id"], g["accession"])
        if not c:
            missing_coords.append(g["gene_id"])
            unplaced.append({**g, "coords": None})
            continue
        sg = (c.get("subgenome") or "").upper()
        row = {**g, "coords": c}
        if sg in by_sg:
            by_sg[sg].append(row)
        else:
            unplaced.append(row)

    rows_out = []
    for sg in ("A", "B", "D"):
        rows_out.append({"subgenome": sg, "entries": by_sg[sg]})
    if unplaced:
        rows_out.append({"subgenome": "Other / no coords", "entries": unplaced})

    return {
        "subgenome_rows": rows_out,
        "missing_coords": sorted(set(missing_coords)),
        "total_in_og": total,
        "shown": len(og_genes),
        "truncated": total > limit,
    }


def build_wheat_homeologue_table(conn, hog_id: str) -> dict:
    """A/B/D rows for genes in the same pangene cluster (requires coords for subgenome)."""
    cur = conn.cursor()
    cur.execute(
        """
        SELECT gene_id, hog, og, accession
        FROM genes WHERE hog = ? ORDER BY accession, gene_id
        """,
        (hog_id,),
    )
    hog_genes = [dict(r) for r in cur.fetchall()]
    by_sg: dict[str, list[dict]] = {"A": [], "B": [], "D": []}
    missing_coords: list[str] = []
    unplaced: list[dict] = []

    for g in hog_genes:
        c = fetch_coords_row(cur, g["gene_id"], g["accession"])
        if not c:
            missing_coords.append(g["gene_id"])
            unplaced.append({**g, "coords": None})
            continue
        sg = (c.get("subgenome") or "").upper()
        row = {**g, "coords": c}
        if sg in by_sg:
            by_sg[sg].append(row)
        else:
            unplaced.append(row)

    rows_out = []
    for sg in ("A", "B", "D"):
        rows_out.append({"subgenome": sg, "entries": by_sg[sg]})
    if unplaced:
        rows_out.append({"subgenome": "Other / no coords", "entries": unplaced})

    return {
        "subgenome_rows": rows_out,
        "missing_coords": sorted(set(missing_coords)),
    }


def build_wheat_og_pangene_triad_panel(
    conn,
    og_id: str,
    current_pangene_id: str,
    *,
    max_hogs: int = 40,
) -> dict | None:
    """
    List distinct pangenes in the same OG, with dominant subgenome (from gene_coords).
    When there are exactly 3 pangenes whose dominants are A, B, D, set triad_columns for a 3-column layout.
    """
    if not table_exists(conn, "gene_coords"):
        return None
    cur = conn.cursor()
    cur.execute("SELECT COUNT(DISTINCT hog) AS n FROM genes WHERE og = ?", (og_id,))
    nh = int(cur.fetchone()["n"])
    if nh > max_hogs:
        return {
            "og": og_id,
            "skipped": True,
            "reason": (
                f"This orthogroup has {nh} distinct pangenes (showing this panel only when ≤ {max_hogs})."
            ),
            "pangene_count": nh,
        }

    cur.execute(
        "SELECT DISTINCT hog FROM genes WHERE og = ? ORDER BY hog",
        (og_id,),
    )
    hogs = [r["hog"] for r in cur.fetchall()]
    rows: list[dict] = []
    for h in hogs:
        cur.execute(
            "SELECT gene_id, accession FROM genes WHERE hog = ?",
            (h,),
        )
        gene_rows = cur.fetchall()
        sub_counts = {"A": 0, "B": 0, "D": 0, "?": 0}
        for g in gene_rows:
            c = fetch_coords_row(cur, g["gene_id"], g["accession"])
            if c and c.get("subgenome"):
                sg = str(c["subgenome"]).upper()
                if sg in sub_counts:
                    sub_counts[sg] += 1
                else:
                    sub_counts["?"] += 1
            else:
                sub_counts["?"] += 1
        best, bestn = "?", -1
        for sg in ("A", "B", "D"):
            if sub_counts[sg] > bestn:
                bestn = sub_counts[sg]
                best = sg
        if bestn == 0:
            best = "?"
        rows.append(
            {
                "pangene": h,
                "dominant": best,
                "counts": sub_counts,
                "gene_count": len(gene_rows),
                "is_current": h == current_pangene_id,
            }
        )

    rows.sort(
        key=lambda r: (
            {"A": 0, "B": 1, "D": 2, "?": 3}.get(r["dominant"], 4),
            r["pangene"],
        )
    )

    triad_columns = None
    if len(rows) == 3:
        doms = sorted(r["dominant"] for r in rows)
        if doms == ["A", "B", "D"]:
            triad_columns = {r["dominant"]: r for r in rows}

    return {
        "og": og_id,
        "skipped": False,
        "rows": rows,
        "triad_columns": triad_columns,
    }


def build_wheat_synteny(
    conn,
    focal_gene_id: str,
    focal_accession: str,
    *,
    color_by: str = "hog",
    window: int = 5,
) -> dict | None:
    """
    Up to (2*window + 1) neighboring genes on same chr/accession (genomic order only).
    SVG uses equal-width arrows (order cartoon, not genomic scale). color by pangene (pan) or OG/type.
    Default window=5 → 11 genes.
    """
    color_by = (color_by or "hog").lower()
    if color_by not in ("hog", "og"):
        color_by = "hog"
    window = max(0, min(int(window), 25))

    cur = conn.cursor()
    fc = fetch_coords_row(cur, focal_gene_id, focal_accession)
    if not fc:
        return None

    chrom = fc["chr"]
    acc = fc["accession"]
    focal_keys = {_norm_gid(x) for x in gene_id_coord_variants(focal_gene_id)}

    slot_count = 2 * window + 1
    slot_js: list[int | None] = []
    chrom_genes: list[dict] = []

    gcm = _gene_coords_meta(conn)
    ci = fc.get("chrom_index") if "chrom_index" in fc else None
    use_window = gcm["window_ok"] and ci is not None and str(ci).strip() != ""
    if use_window:
        try:
            F = int(ci)
        except (TypeError, ValueError):
            use_window = False

    if use_window:
        cur.execute(
            """
            SELECT MAX(chrom_index) AS m FROM gene_coords
            WHERE accession = ? AND chr = ?
            """,
            (acc, chrom),
        )
        mr = cur.fetchone()
        max_ci = int(mr["m"]) if mr and mr["m"] is not None else F
        first_ci = max(0, F - window)
        last_ci = min(max_ci, F + window)
        cur.execute(
            """
            SELECT gene_id, chr, start, end, strand, subgenome, chrom_index
            FROM gene_coords
            WHERE accession = ? AND chr = ?
              AND chrom_index BETWEEN ? AND ?
            ORDER BY chrom_index ASC
            """,
            (acc, chrom, first_ci, last_ci),
        )
        by_ci = {int(r["chrom_index"]): dict(r) for r in cur.fetchall()}
        # Verify focal is at F (variant IDs / stale index guard).
        g_at_f = by_ci.get(F)
        if not g_at_f or _norm_gid(g_at_f["gene_id"]) not in focal_keys:
            use_window = False
            slot_js = []
            chrom_genes = []

    if use_window:
        for k in range(slot_count):
            tci = F - window + k
            if tci < 0 or tci > max_ci:
                slot_js.append(None)
            elif tci in by_ci:
                slot_js.append(len(chrom_genes))
                chrom_genes.append(by_ci[tci])
            else:
                slot_js.append(None)
    else:
        cur.execute(
            """
            SELECT gene_id, chr, start, end, strand, subgenome
            FROM gene_coords
            WHERE accession = ? AND chr = ?
            ORDER BY start ASC, end ASC
            """,
            (acc, chrom),
        )
        chrom_genes = [dict(r) for r in cur.fetchall()]
        if not chrom_genes:
            return None

        idx = None
        for i, g in enumerate(chrom_genes):
            if _norm_gid(g["gene_id"]) in focal_keys:
                idx = i
                break
        if idx is None:
            return None

        for k in range(slot_count):
            j = idx - window + k
            slot_js.append(j if 0 <= j < len(chrom_genes) else None)

    if not chrom_genes and not any(j is not None for j in slot_js):
        return None

    present = [chrom_genes[j] for j in slot_js if j is not None]
    if not present:
        return None
    min_s = min(g["start"] for g in present)
    max_e = max(g["end"] for g in present)

    ortho_batch = fetch_gene_ortho_batch(
        cur, [g["gene_id"] for g in present]
    )

    missing_ortho: list[str] = []
    segments: list[dict] = []
    legend = []
    seen_leg: set[str] = set()

    for j in slot_js:
        if j is None:
            segments.append({"empty": True})
            continue
        g = chrom_genes[j]
        ortho = ortho_batch.get(g["gene_id"])
        if not ortho:
            missing_ortho.append(g["gene_id"])
            ck = "—"
        else:
            ck = ortho["hog"] if color_by == "hog" else ortho["og"]
        is_focal = _norm_gid(g["gene_id"]) in focal_keys
        strand = (g.get("strand") or "").lower()
        color = _synteny_color_for_key(ck)
        segments.append(
            {
                "empty": False,
                "strand": strand,
                "color": color,
                "gene_id": g["gene_id"],
                "accession": acc,
                "chr": g["chr"],
                "start": int(g["start"]),
                "end": int(g["end"]),
                "is_focal": is_focal,
                "pangene": ortho["hog"] if ortho else "",
                "og": ortho["og"] if ortho else "",
                "color_key": ck,
            }
        )
        if ck not in seen_leg:
            seen_leg.add(ck)
            legend.append({"key": ck, "color": color})

    svg_w = 900
    svg_h = 108
    return {
        "accession": acc,
        "chrom": chrom,
        "color_by": color_by,
        "min_pos": min_s,
        "max_pos": max_e,
        "svg_w": svg_w,
        "svg_h": svg_h,
        "segments": segments,
        "legend": legend,
        "missing_ortho": sorted(set(missing_ortho)),
        "equal_width": True,
    }


def render_wheat_synteny_svg(syn: dict) -> str:
    """Equal-width arrow row: genomic order only, not gene length or intergenic distance."""
    segs = syn.get("segments") or []
    if not segs:
        return ""
    compact = bool(syn.get("compact"))
    w = float(syn["svg_w"])
    n = len(segs)
    margin = 8.0 if compact else 16.0
    gap = 4.0 if compact else 5.0
    if compact:
        h = 48.0
        body_h = 22.0
        y_c = h / 2
        tip_frac = 0.22
    else:
        h = float(syn["svg_h"])
        body_h = 28.0
        y_c = 58.0
        tip_frac = 0.2
    inner = w - 2 * margin
    total_gap = gap * max(0, n - 1)
    slot_w = max((inner - total_gap) / n, 8.0 if compact else 6.0)
    tip_w = min(max(slot_w * tip_frac, 7.0 if compact else 6.0), slot_w * 0.45)
    body_w = max(slot_w - tip_w, 4.0)
    y0 = y_c - body_h / 2
    y1 = y_c + body_h / 2

    cls = "synteny-svg synteny-svg-equal synteny-svg-compact" if compact else "synteny-svg synteny-svg-equal"
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" '
        f'class="{cls}" role="img" aria-label="Gene order along chromosome">'
    ]
    if not compact:
        n_real = sum(1 for x in segs if not x.get("empty"))
        sub = (
            f"{escape(syn['accession'])} chr{escape(str(syn['chrom']))} · "
            f"{n_real} genes in {n} slots · equal-width (order only) · color "
            f"{escape(syn['color_by'].upper())}"
        )
        parts.append(
            f'<text x="12" y="20" font-size="11" fill="#555">{sub}</text>'
        )
        parts.append(
            f'<text x="12" y="36" font-size="10" fill="#888">'
            f"Locus ~{syn['min_pos']:,}–{syn['max_pos']:,} bp (not to scale here)</text>"
        )

    for i, s in enumerate(segs):
        x_slot = margin + i * (slot_w + gap)
        if s.get("empty"):
            parts.append(
                f'<rect x="{x_slot:.1f}" y="{y0:.1f}" width="{slot_w:.1f}" '
                f'height="{body_h:.1f}" rx="2" fill="#e9ecef" stroke="#ced4da" '
                f'stroke-width="0.8" pointer-events="none"/>'
            )
            continue
        strand = (s.get("strand") or "").lower()
        col = s["color"]
        stroke = "#1a1a1a" if s.get("is_focal") else "#555"
        sw = 2.2 if s.get("is_focal") else 1.1
        acc_raw = str(s.get("accession") or syn.get("accession") or "")
        loc_plain = (
            f"{s.get('chr') or ''}:{int(s.get('start') or 0):,}-{int(s.get('end') or 0):,}"
        )
        parts.append(
            '<g class="synteny-seg" '
            f'data-accession="{escape(acc_raw)}" '
            f'data-gene-id="{escape(s["gene_id"])}" '
            f'data-locus="{escape(loc_plain)}" '
            f'data-pangene="{escape(str(s.get("pangene") or ""))}" '
            f'data-og="{escape(str(s.get("og") or ""))}">'
        )

        if strand in ("-", "rev", "reverse", "r"):
            x_tip = x_slot
            x_rect = x_tip + tip_w
            pts = (
                f"{x_tip:.1f},{y_c:.1f} {x_rect:.1f},{y0:.1f} "
                f"{x_rect + body_w:.1f},{y0:.1f} {x_rect + body_w:.1f},{y1:.1f} "
                f"{x_rect:.1f},{y1:.1f}"
            )
        else:
            x_rect = x_slot
            x_tip = x_rect + body_w + tip_w
            pts = (
                f"{x_rect:.1f},{y0:.1f} {x_rect + body_w:.1f},{y0:.1f} "
                f"{x_tip:.1f},{y_c:.1f} {x_rect + body_w:.1f},{y1:.1f} "
                f"{x_rect:.1f},{y1:.1f}"
            )

        parts.append(
            f'<polygon points="{pts}" fill="{col}" stroke="{stroke}" stroke-width="{sw}" '
            f'stroke-linejoin="round" pointer-events="visiblePainted"/>'
        )
        parts.append("</g>")

    parts.append("</svg>")
    return "".join(parts)


def enrich_wheat_results_groups(hog_groups: list[dict], query: str, color_by: str):
    conn = get_db("wheat")
    try:
        if not table_exists(conn, "gene_coords"):
            for g in hog_groups:
                g["wheat_homeologue"] = None
                g["wheat_og_subgenome"] = None
                g["wheat_og_pangene_triad"] = None
                g["wheat_synteny"] = None
                g["wheat_synteny_svg"] = ""
            return
        qn = query.strip().lower()

        for group in hog_groups:
            pan = group["pangene"]
            group["wheat_homeologue"] = build_wheat_homeologue_table(conn, pan)
            group["wheat_og_subgenome"] = build_wheat_og_subgenome_table(
                conn, group["og"]
            )
            group["wheat_og_pangene_triad"] = build_wheat_og_pangene_triad_panel(
                conn, group["og"], pan
            )

            focal = None
            for ge in group["genes"]:
                if ge["gene_id"].lower() == qn:
                    focal = ge
                    break
            if focal is None:
                for ge in group["genes"]:
                    if qn in ge["gene_id"].lower():
                        focal = ge
                        break
            if focal is None and group["genes"]:
                focal = group["genes"][0]

            syn = None
            if focal:
                syn = build_wheat_synteny(
                    conn,
                    focal["gene_id"],
                    focal["accession"],
                    color_by=color_by,
                    window=5,
                )
            group["wheat_synteny"] = syn
            group["wheat_synteny_svg"] = render_wheat_synteny_svg(syn) if syn else ""
    finally:
        conn.close()


def enrich_wheat_pandagma_results_groups(
    hog_groups: list[dict], query: str, color_by: str
):
    """
    Pandagma mode enrichment:
    - keep only local synteny for the query match window
    - avoid OG/pangene-triad panels that assume OrthoFinder-style sibling clusters
    """
    conn = get_db("wheat")
    try:
        if not table_exists(conn, "gene_coords"):
            for g in hog_groups:
                g["wheat_homeologue"] = None
                g["wheat_og_subgenome"] = None
                g["wheat_og_pangene_triad"] = None
                g["wheat_synteny"] = None
                g["wheat_synteny_svg"] = ""
            return

        qn = query.strip().lower()
        for group in hog_groups:
            group["wheat_homeologue"] = None
            group["wheat_og_subgenome"] = None
            group["wheat_og_pangene_triad"] = None

            focal = None
            for ge in group["genes"]:
                if ge["gene_id"].lower() == qn:
                    focal = ge
                    break
            if focal is None:
                for ge in group["genes"]:
                    if qn in ge["gene_id"].lower():
                        focal = ge
                        break
            if focal is None and group["genes"]:
                focal = group["genes"][0]

            syn = None
            if focal:
                syn = build_wheat_synteny(
                    conn,
                    focal["gene_id"],
                    focal["accession"],
                    color_by=color_by,
                    window=5,
                )
            group["wheat_synteny"] = syn
            group["wheat_synteny_svg"] = render_wheat_synteny_svg(syn) if syn else ""
    finally:
        conn.close()


# --- Routes ---

# Supported datasets are those with both a DB and FASTA dir present.
DATASETS_WITH_DATA = available_datasets()


def default_dataset_id() -> str:
    """Prefer wheat when available; otherwise fall back to any installed dataset."""
    if "wheat" in DATASETS_WITH_DATA:
        return "wheat"
    if DATASETS_WITH_DATA:
        return sorted(DATASETS_WITH_DATA)[0]
    return DEFAULT_DATASET_ID


@app.context_processor
def inject_nav_defaults():
    return {"nav_default_dataset": default_dataset_id()}


@app.route("/dataset/<dataset_id>")
def dataset(dataset_id):
    dataset_id = (dataset_id or "").strip().lower()
    if dataset_id in DATASETS_WITH_DATA:
        from flask import redirect, url_for
        return redirect(url_for("index") + f"?dataset_id={dataset_id}")
    return render_template(
        "dataset_unavailable.html",
        dataset_id=dataset_id,
        dataset_label=dataset_id.capitalize(),
    )


@app.route("/")
def index():
    dataset_id = request.args.get("dataset_id", default_dataset_id()).strip().lower()
    if dataset_id not in DATASET_CONFIG:
        dataset_id = default_dataset_id()
    return render_template(
        "index.html",
        dataset_id=dataset_id,
        dataset_label=DATASET_CONFIG[dataset_id]["label"],
        about_stats=load_about_stats(dataset_id),
    )


@app.route("/search")
def search():
    dataset_id = request.args.get("dataset_id", default_dataset_id()).strip().lower()
    if dataset_id not in DATASET_CONFIG:
        dataset_id = default_dataset_id()
    query = request.args.get("q", "").strip()
    if not query:
        return render_template(
            "index.html",
            error="Please enter a gene ID.",
            dataset_id=dataset_id,
            dataset_label=DATASET_CONFIG[dataset_id]["label"],
            about_stats=load_about_stats(dataset_id),
        )

    results = search_genes(query, dataset_id)
    if not results:
        return render_template(
            "index.html",
            error=f'No results found for "{query}".',
            query=query,
            dataset_id=dataset_id,
            dataset_label=DATASET_CONFIG[dataset_id]["label"],
            about_stats=load_about_stats(dataset_id),
        )

    if len(results) == 1:
        gene = results[0]
        from flask import redirect, url_for

        return redirect(
            url_for(
                "pangene_detail",
                pangene_id=gene["hog"],
                highlight=query,
                dataset_id=dataset_id,
            )
        )

    pans_seen: dict[str, dict] = {}
    for r in results:
        pan = r["hog"]
        if pan not in pans_seen:
            pans_seen[pan] = {
                "pangene": pan,
                "og": r["og"],
                "genes": [],
            }
        pans_seen[pan]["genes"].append(
            {"gene_id": r["gene_id"], "accession": r["accession"]}
        )

    pangene_list = list(pans_seen.values())
    synteny_color = request.args.get("synteny_color", "hog").strip().lower()
    if synteny_color not in ("hog", "og"):
        synteny_color = "hog"

    wheat_coords_available = False
    if dataset_id == "wheat":
        _wc = get_db("wheat")
        wheat_coords_available = table_exists(_wc, "gene_coords")
        _wc.close()
        if wheat_coords_available:
            cfg = get_dataset_config(dataset_id)
            if cfg.get("mode") == "pandagma" and not cfg.get(
                "enable_homeologue_panels"
            ):
                enrich_wheat_pandagma_results_groups(pangene_list, query, synteny_color)
            else:
                enrich_wheat_results_groups(pangene_list, query, synteny_color)

    return render_template(
        "results.html",
        query=query,
        dataset_id=dataset_id,
        pangene_groups=pangene_list,
        synteny_color=synteny_color,
        show_wheat_extra=(dataset_id == "wheat"),
        wheat_coords_available=wheat_coords_available,
    )


@app.route("/pangene/<pangene_id>")
def pangene_detail(pangene_id):
    dataset_id = request.args.get("dataset_id", default_dataset_id()).strip().lower()
    if dataset_id not in DATASET_CONFIG:
        dataset_id = default_dataset_id()
    highlight = request.args.get("highlight", "")
    info = get_hog_info(pangene_id, dataset_id)
    if not info:
        return render_template(
            "index.html",
            error=f'Pangene "{pangene_id}" not found.',
            dataset_id=dataset_id,
            dataset_label=DATASET_CONFIG[dataset_id]["label"],
            about_stats=load_about_stats(dataset_id),
        )

    genes = get_hog_genes(pangene_id, dataset_id)
    sequences = read_fasta(pangene_id, dataset_id)

    aligned = {}
    consensus_str = ""
    aln_len = 0

    newick = ""
    leaf_order = []

    conservation_scores = []
    if sequences and len(sequences) > 0:
        aligned, newick = run_famsa(sequences)
        if aligned:
            consensus_str = compute_consensus(aligned)
            aln_len = max(len(s) for s in aligned.values())
            leaf_order = get_tree_leaf_order(newick)
            conservation_scores = compute_conservation(aligned)

    gene_accession_lookup = build_gene_accession_lookup(pangene_id, dataset_id)

    wheat_og_picker = None
    wheat_synteny_available = False
    gene_row_meta: dict[str, dict] = {}
    pandagma_pan_view_available = False
    cds_available = False
    porter6_by_gene: dict[str, dict[str, str]] = {}
    _conn = get_db(dataset_id)
    try:
        cds_available = table_exists(_conn, "cds_seqs")
        wheat_synteny_available = table_exists(_conn, "gene_coords")
        cfg = get_dataset_config(dataset_id)
        pandagma_pan_view_available = (
            dataset_id == "wheat"
            and cfg.get("mode") == "pandagma"
            and wheat_synteny_available
        )
        if dataset_id == "wheat" and get_dataset_config(dataset_id).get(
            "enable_og_picker"
        ):
            wheat_og_picker = build_wheat_og_picker_payload(
                _conn, info["og"], pangene_id
            )
        cur = _conn.cursor()
        coords_on = wheat_synteny_available
        ga_pairs = [(row["gene_id"], row["accession"]) for row in genes]
        coords_by_gid = (
            fetch_coords_rows_batch(cur, ga_pairs) if coords_on else {}
        )
        for row in genes:
            gid = row["gene_id"]
            acc = row["accession"]
            m = {
                "og": info["og"],
                "pangene": pangene_id,
                "chr": None,
                "start": None,
                "end": None,
                "strand": None,
                "subgenome": None,
            }
            if coords_on:
                cr = coords_by_gid.get(gid)
                if cr:
                    m["chr"] = cr["chr"]
                    m["start"] = int(cr["start"])
                    m["end"] = int(cr["end"])
                    st = (cr.get("strand") or "").strip()
                    m["strand"] = st if st else None
                    sg = (cr.get("subgenome") or "").strip().upper()
                    m["subgenome"] = sg if sg in ("A", "B", "D") else None
            gene_row_meta[gid] = m
        porter6_by_gene = (
            porter6_raw_for_aligned_sequences(_conn, aligned) if aligned else {}
        )
    finally:
        _conn.close()

    chinese_spring_gene_ids = chinese_spring_gene_ids_from_rows(genes)

    return render_template(
        "pangene_detail.html",
        dataset_id=dataset_id,
        pangene_id=pangene_id,
        info=info,
        genes=genes,
        gene_row_meta=gene_row_meta,
        sequences=sequences,
        aligned_json=json.dumps(aligned),
        consensus_str=consensus_str,
        aln_len=aln_len,
        highlight=highlight,
        num_seqs=len(aligned),
        newick=newick,
        leaf_order_json=json.dumps(leaf_order),
        conservation_json=json.dumps(conservation_scores),
        wheat_og_picker=wheat_og_picker,
        wheat_synteny_available=wheat_synteny_available,
        gene_accession_lookup=gene_accession_lookup,
        pandagma_pan_view_available=pandagma_pan_view_available,
        cds_available=cds_available,
        porter6_json=porter6_by_gene,
        external_db_links=EXTERNAL_DB_LINKS,
        chinese_spring_gene_ids=chinese_spring_gene_ids,
    )


@app.route("/api/plantapp_omics")
def api_plantapp_omics():
    """
    Proxy PlantApp tissue + DEG JSON for one gene_id (see PlantApp pages/api.py).
    """
    gene_id = (request.args.get("gene_id") or "").strip()
    if not gene_id:
        return jsonify({"error": "gene_id is required"}), 400
    return jsonify(fetch_plantapp_omics(gene_id))


@app.route("/pangene/<pangene_id>/tree")
def pangene_tree(pangene_id):
    dataset_id = request.args.get("dataset_id", default_dataset_id()).strip().lower()
    if dataset_id not in DATASET_CONFIG:
        dataset_id = default_dataset_id()
    info = get_hog_info(pangene_id, dataset_id)
    if not info:
        return render_template(
            "index.html",
            error=f'Pangene "{pangene_id}" not found.',
            dataset_id=dataset_id,
            dataset_label=DATASET_CONFIG[dataset_id]["label"],
            about_stats=load_about_stats(dataset_id),
        )

    sequences = read_fasta(pangene_id, dataset_id)
    aligned = {}
    consensus_str = ""
    newick = ""
    aln_len = 0
    leaf_order = []

    if sequences and len(sequences) > 1:
        aligned, newick = run_famsa(sequences)
        if aligned:
            consensus_str = compute_consensus(aligned)
            aln_len = max(len(s) for s in aligned.values())
            leaf_order = get_tree_leaf_order(newick)

    return render_template(
        "tree.html",
        dataset_id=dataset_id,
        pangene_id=pangene_id,
        info=info,
        newick=newick,
        aligned_json=json.dumps(aligned),
        consensus_str=consensus_str,
        aln_len=aln_len,
        leaf_order_json=json.dumps(leaf_order),
        num_seqs=len(aligned),
    )


def compute_variants_for_download(
    aligned: dict[str, str],
    *,
    coords_by_gid: dict[str, dict | None],
    acc_by_gid: dict[str, str],
) -> str:
    """
    Transposed TSV: columns are alignment positions (1-based) where any gene differs
    from the first gene (alphabetically); rows are genes with gene_id, accession, chr,
    start, end, then one AA per variant column.
    """
    if not aligned:
        return ""
    gids = sorted(aligned.keys())
    aln_len = max(len(aligned[g]) for g in gids)
    ref = gids[0]
    ref_seq = aligned[ref]
    variant_positions: list[int] = []
    for i in range(aln_len):
        ra = ref_seq[i] if i < len(ref_seq) else "-"
        for g in gids:
            aa = aligned[g][i] if i < len(aligned[g]) else "-"
            if aa != ra:
                variant_positions.append(i)
                break
    headers = ["gene_id", "accession", "chr", "start", "end"] + [
        str(p + 1) for p in variant_positions
    ]
    lines = ["\t".join(headers)]
    for gid in gids:
        acc = acc_by_gid.get(gid, "") or ""
        cr = coords_by_gid.get(gid)
        chr_s = ""
        start_s = ""
        end_s = ""
        if cr:
            if cr.get("chr") is not None:
                chr_s = str(cr["chr"])
            if cr.get("start") is not None:
                start_s = str(int(cr["start"]))
            if cr.get("end") is not None:
                end_s = str(int(cr["end"]))
        row_cells = [gid, acc, chr_s, start_s, end_s]
        for p in variant_positions:
            aa = aligned[gid][p] if p < len(aligned[gid]) else "-"
            row_cells.append(aa if aa not in (".",) else "-")
        lines.append("\t".join(row_cells))
    return "\n".join(lines)


@app.route("/download/<pangene_id>/<dtype>")
def download(pangene_id, dtype):
    dataset_id = request.args.get("dataset_id", default_dataset_id()).strip().lower()
    if dataset_id not in DATASET_CONFIG:
        dataset_id = default_dataset_id()

    info = get_hog_info(pangene_id, dataset_id)
    if not info:
        return "Pangene not found", 404

    if dtype == "gene_list":
        genes = get_hog_genes(pangene_id, dataset_id)
        lines = ["gene_id\taccession\tchromosome\tstart\tend"]
        for g in genes:
            row = dict(g)
            gid = row["gene_id"]
            acc = row["accession"]
            ch = row.get("chr")
            st = row.get("start")
            en = row.get("end")
            lines.append(
                "\t".join(
                    [
                        gid,
                        acc,
                        "" if ch is None else str(ch),
                        "" if st is None else str(int(st)),
                        "" if en is None else str(int(en)),
                    ]
                )
            )
        content = "\n".join(lines)
        return Response(
            content,
            mimetype="text/plain",
            headers={
                "Content-Disposition": f"attachment; filename={pangene_id}_genes.tsv"
            },
        )

    elif dtype == "sequences":
        sequences = read_fasta(pangene_id, dataset_id)
        lines = []
        for sid, seq in sequences.items():
            lines.append(f">{sid}")
            for i in range(0, len(seq), 80):
                lines.append(seq[i : i + 80])
        content = "\n".join(lines)
        return Response(
            content,
            mimetype="text/plain",
            headers={
                "Content-Disposition": f"attachment; filename={pangene_id}_proteins.fasta"
            },
        )

    elif dtype == "alignment":
        sequences = read_fasta(pangene_id, dataset_id)
        if sequences:
            aligned, _ = run_famsa(sequences)
            lines = []
            for sid, seq in aligned.items():
                lines.append(f">{sid}")
                for i in range(0, len(seq), 80):
                    lines.append(seq[i : i + 80])
            content = "\n".join(lines)
        else:
            content = ""
        return Response(
            content,
            mimetype="text/plain",
            headers={
                "Content-Disposition": f"attachment; filename={pangene_id}_alignment.fasta"
            },
        )

    elif dtype == "variants":
        sequences = read_fasta(pangene_id, dataset_id)
        if sequences:
            aligned, _ = run_famsa(sequences)
            acc_by_gid = build_gene_accession_lookup(pangene_id, dataset_id)
            conn = get_db(dataset_id)
            try:
                coords_on = table_exists(conn, "gene_coords")
                cur = conn.cursor()
                pairs = [(g, acc_by_gid.get(g, "")) for g in sorted(aligned.keys())]
                coords_by_gid = (
                    fetch_coords_rows_batch(cur, pairs) if coords_on else {}
                )
            finally:
                conn.close()
            content = compute_variants_for_download(
                aligned,
                coords_by_gid=coords_by_gid,
                acc_by_gid=acc_by_gid,
            )
        else:
            content = ""
        return Response(
            content,
            mimetype="text/plain",
            headers={
                "Content-Disposition": f"attachment; filename={pangene_id}_variants.tsv"
            },
        )

    return "Unknown download type", 400


@app.route("/api/search_suggestions")
def search_suggestions():
    dataset_id = request.args.get("dataset_id", default_dataset_id()).strip().lower()
    if dataset_id not in DATASET_CONFIG:
        dataset_id = default_dataset_id()
    query = request.args.get("q", "").strip()
    if len(query) < 3:
        return jsonify([])
    conn = get_db(dataset_id)
    cur = conn.cursor()
    cur.execute(
        "SELECT gene_id, hog, accession FROM genes "
        "WHERE gene_id LIKE ? COLLATE NOCASE LIMIT 15",
        (f"%{query}%",),
    )
    results = [
        {
            "gene_id": r["gene_id"],
            "pangene": r["hog"],
            "accession": r["accession"],
        }
        for r in cur.fetchall()
    ]
    conn.close()
    return jsonify(results)


@app.route("/api/wheat/synteny_tracks", methods=["POST"])
def api_wheat_synteny_tracks():
    """One ±window synteny track per pangene (focal = first gene with coords in that cluster)."""
    if "wheat" not in DATASETS_WITH_DATA:
        return jsonify(error="wheat dataset not available"), 404
    conn = get_db("wheat")
    try:
        if not table_exists(conn, "gene_coords"):
            return jsonify(error="gene_coords not loaded"), 503
        data = request.get_json(silent=True) or {}
        raw_ids = data.get("pangene_ids") or []
        color_by = (data.get("color_by") or "hog").strip().lower()
        if color_by not in ("hog", "og"):
            color_by = "hog"
        if not isinstance(raw_ids, list):
            return jsonify(error="pangene_ids must be a list"), 400
        pangene_ids = [str(h).strip() for h in raw_ids if str(h).strip()]
        cur = conn.cursor()
        tracks: list[dict] = []
        for hid in pangene_ids:
            cur.execute(
                "SELECT gene_id, accession FROM genes WHERE hog = ? ORDER BY gene_id",
                (hid,),
            )
            genes = cur.fetchall()
            if not genes:
                tracks.append(
                    {"pangene": hid, "ok": False, "message": "Pangene not found"}
                )
                continue
            focal = None
            for g in genes:
                if fetch_coords_row(cur, g["gene_id"], g["accession"]):
                    focal = dict(g)
                    break
            if not focal:
                tracks.append(
                    {
                        "pangene": hid,
                        "ok": False,
                        "message": "No gene with coordinates in gene_coords",
                    }
                )
                continue
            syn = build_wheat_synteny(
                conn,
                focal["gene_id"],
                focal["accession"],
                color_by=color_by,
            )
            if not syn:
                tracks.append(
                    {
                        "pangene": hid,
                        "ok": False,
                        "focal_gene_id": focal["gene_id"],
                        "message": "Could not build synteny window",
                    }
                )
                continue
            tracks.append(
                {
                    "pangene": hid,
                    "ok": True,
                    "focal_gene_id": focal["gene_id"],
                    "focal_accession": focal["accession"],
                    "svg": render_wheat_synteny_svg(syn),
                    "legend": syn.get("legend", []),
                    "missing_ortho": syn.get("missing_ortho", []),
                }
            )
        return jsonify(tracks=tracks)
    finally:
        conn.close()


@app.route("/api/wheat/merged_alignment", methods=["POST"])
def api_wheat_merged_alignment():
    """FAMSA (NJ guide tree) across proteins for the listed pangene clusters (wheat)."""
    if "wheat" not in DATASETS_WITH_DATA:
        return jsonify(error="wheat dataset not available"), 404
    data = request.get_json(silent=True) or {}
    raw = data.get("pangene_ids") or []
    pangene_ids = sorted({str(h).strip() for h in raw if str(h).strip()})
    if len(pangene_ids) < 1:
        return jsonify(error="no_pangene_ids"), 400

    merged: dict[str, str] = {}
    missing_fastas: list[str] = []
    for hid in pangene_ids:
        seqs = read_fasta(hid, "wheat")
        if not seqs:
            missing_fastas.append(hid)
        merged.update(seqs)

    if len(merged) < 2:
        return jsonify(
            error="not_enough_sequences",
            missing_fastas=missing_fastas,
        ), 400

    aligned, newick = run_famsa(merged)
    aln_len = max((len(s) for s in aligned.values()), default=0)
    lo = get_tree_leaf_order(newick)
    cons = compute_conservation(aligned)
    acc_map = build_multi_hog_gene_accession_lookup(pangene_ids, "wheat")

    conn = get_db("wheat")
    try:
        porter6 = porter6_raw_for_aligned_sequences(conn, aligned)
    finally:
        conn.close()

    return jsonify(
        aligned=aligned,
        aln_len=aln_len,
        newick=newick or "",
        leaf_order=lo,
        conservation=cons,
        gene_accession=acc_map,
        porter6=porter6,
        num_seqs=len(aligned),
        missing_fastas=missing_fastas,
    )


@app.route("/api/align_genes", methods=["POST"])
def api_align_genes():
    """Re-align only the requested gene IDs (active subset) with FAMSA."""
    data = request.get_json(silent=True) or {}
    dataset_id = (data.get("dataset_id") or default_dataset_id()).strip().lower()
    if dataset_id not in DATASET_CONFIG:
        return jsonify(error="unknown_dataset"), 400
    if dataset_id not in DATASETS_WITH_DATA:
        return jsonify(error="dataset_not_available"), 404
    raw_ids = data.get("gene_ids") or []
    if not isinstance(raw_ids, list):
        return jsonify(error="gene_ids_must_be_list"), 400
    seen: set[str] = set()
    gene_ids: list[str] = []
    for g in raw_ids:
        s = str(g).strip()
        if s and s not in seen:
            seen.add(s)
            gene_ids.append(s)
    if len(gene_ids) < 1:
        return jsonify(error="no_genes"), 400
    if len(gene_ids) > 600:
        return jsonify(error="too_many_genes_max_600"), 400

    conn = get_db(dataset_id)
    try:
        merged, acc_map = collect_sequences_for_gene_ids(conn, gene_ids, dataset_id)
    finally:
        conn.close()

    if len(merged) < 1:
        return jsonify(error="no_sequences_found"), 400

    if len(merged) == 1:
        only_g, only_s = next(iter(merged.items()))
        aligned = {only_g: only_s}
        aln_len = len(only_s)
        newick = ""
        lo = [only_g]
        cons = [1.0] * aln_len
    else:
        aligned, newick = run_famsa(merged)
        aln_len = max((len(s) for s in aligned.values()), default=0)
        lo = get_tree_leaf_order(newick)
        cons = compute_conservation(aligned)

    out_acc = {g: acc_map.get(g, g) for g in aligned}

    conn2 = get_db(dataset_id)
    try:
        porter6 = porter6_raw_for_aligned_sequences(conn2, aligned)
    finally:
        conn2.close()

    return jsonify(
        aligned=aligned,
        aln_len=aln_len,
        newick=newick or "",
        leaf_order=lo,
        conservation=cons,
        gene_accession=out_acc,
        porter6=porter6,
        num_seqs=len(aligned),
    )


@app.route("/api/ss_plot", methods=["POST"])
def api_ss_plot():
    """Return an SVG secondary-structure diagram (secstructartist + matplotlib; no temp files)."""
    data = request.get_json(silent=True) or {}
    raw_rows = data.get("rows")
    mode = (data.get("mode") or "q3").strip().lower()
    if not isinstance(raw_rows, list):
        return jsonify(error="rows_must_be_list"), 400
    if mode not in ("q3", "q8"):
        return jsonify(error="bad_mode"), 400

    clean: list[str] = []
    total_len = 0
    for r in raw_rows[:120]:
        if not isinstance(r, str):
            continue
        s = r[:50]
        if not s or not any(ch != " " for ch in s):
            continue
        clean.append(s)
        total_len += len(s)
    if not clean:
        return jsonify(error="no_rows"), 400
    if total_len > 6000:
        return jsonify(error="too_large"), 400

    svg, err = secstructartist_svg_from_rows(clean, mode)
    if svg is None:
        return jsonify(error=err or "render_failed"), 500
    return jsonify(svg=svg)


@app.route("/api/wheat/synteny_per_gene", methods=["POST"])
def api_wheat_synteny_per_gene():
    """One compact equal-width synteny strip per (gene_id, accession) focal."""
    if "wheat" not in DATASETS_WITH_DATA:
        return jsonify(error="wheat dataset not available"), 404
    conn = get_db("wheat")
    try:
        if not table_exists(conn, "gene_coords"):
            return jsonify(error="gene_coords not loaded"), 503
        data = request.get_json(silent=True) or {}
        genes_in = data.get("genes") or []
        color_by = (data.get("color_by") or "hog").strip().lower()
        if color_by not in ("hog", "og"):
            color_by = "hog"
        if not isinstance(genes_in, list):
            return jsonify(error="genes must be a list"), 400
        if len(genes_in) > 150:
            return jsonify(error="too_many_genes_max_150"), 400

        tracks: list[dict] = []
        cur = conn.cursor()
        for item in genes_in:
            if not isinstance(item, dict):
                continue
            gid = str(item.get("gene_id") or "").strip()
            acc = str(item.get("accession") or "").strip()
            if not gid or not acc:
                continue
            fc = fetch_coords_row(cur, gid, acc)
            syn = build_wheat_synteny(
                conn, gid, acc, color_by=color_by, window=5
            )
            if not syn:
                tracks.append(
                    {
                        "gene_id": gid,
                        "accession": acc,
                        "ok": False,
                        "chrom": fc["chr"] if fc else None,
                        "message": "No synteny window (missing coords or locus)",
                        "svg": "",
                    }
                )
                continue
            syn = dict(syn)
            syn["svg_w"] = 720
            syn["compact"] = True
            # Client renders SVG (threshold coloring + session cache without refetch).
            tracks.append(
                {
                    "gene_id": gid,
                    "accession": acc,
                    "ok": True,
                    "chrom": syn["chrom"],
                    "syn": syn,
                    "missing_ortho": syn.get("missing_ortho", []),
                }
            )

        return Response(json.dumps(tracks), mimetype="application/json")
    finally:
        conn.close()


@app.route("/api/pairwise_kaks", methods=["POST"])
def api_pairwise_kaks():
    """
    Pairwise Ka/Ks vs a reference from CDS mapped onto the client-provided protein alignment.
    """
    data = request.get_json(silent=True) or {}
    dataset_id = (data.get("dataset_id") or default_dataset_id()).strip().lower()
    if dataset_id not in DATASET_CONFIG:
        dataset_id = default_dataset_id()
    pangene_id = (data.get("pangene_id") or "").strip()
    ref_gene_id = (data.get("ref_gene_id") or "").strip()
    aligned = data.get("aligned")
    if not pangene_id or not ref_gene_id or not isinstance(aligned, dict):
        return jsonify(error="bad_request"), 400
    if ref_gene_id not in aligned:
        return jsonify(error="ref_not_in_alignment"), 400
    if len(aligned) > 200:
        return jsonify(error="too_many_sequences"), 400
    max_len = max((len(s) for s in aligned.values() if isinstance(s, str)), default=0)
    if max_len > 30000:
        return jsonify(error="alignment_too_long"), 400

    conn = get_db(dataset_id)
    try:
        if not table_exists(conn, "cds_seqs"):
            return jsonify(error="cds_not_loaded"), 503
        cur = conn.cursor()
        cur.execute("SELECT gene_id FROM genes WHERE hog = ?", (pangene_id,))
        allowed = {r["gene_id"] for r in cur.fetchall()}
        for gid in aligned:
            if gid not in allowed:
                return jsonify(error="gene_not_in_group", gene_id=gid), 400
        rows = compute_kaks_vs_reference(conn, pangene_id, ref_gene_id, aligned)
        try:
            window_half = int(data.get("window_half", 5))
        except (TypeError, ValueError):
            window_half = 5
        window_half = max(1, min(window_half, 80))
        profile = compute_sliding_kaks_profile(
            conn, pangene_id, ref_gene_id, aligned, window_half=window_half
        )
        return jsonify(
            ok=True,
            ref_gene_id=ref_gene_id,
            rows=rows,
            profile=profile,
            window_half=window_half,
        )
    finally:
        conn.close()


if __name__ == "__main__":
    _dd = default_dataset_id()
    if not os.path.exists(DATASET_CONFIG[_dd]["db_path"]):
        _dd = DEFAULT_DATASET_ID
    default_db = DATASET_CONFIG[_dd]["db_path"]
    if not os.path.exists(default_db):
        print(
            "Database not found for the default dataset. "
            "Run `python build_index.py` (Pandagma wheat) or "
            "`python build_index.py --n0 <N0.tsv> --out <db_path>` first."
        )
        exit(1)
    app.run(debug=True, port=5050)
