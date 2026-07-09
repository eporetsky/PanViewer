import csv
import hashlib
import io
import json
import logging
import math
import os
import re
import secrets
import time
import shutil
import sqlite3
import subprocess
import tempfile
import threading
from datetime import date
import urllib.error
import urllib.request
import zlib
from collections import OrderedDict, defaultdict
from typing import Any

from markupsafe import escape

from plantapp_omics import fetch_plantapp_omics

from gene_id_normalize import canonical_gene_id, fasta_header_token_candidates
from track_order_kmer import protein_track_order_permutation

from database_config import (
    dataset_configs,
    default_variant_for_species,
    keyword_dataset_id_for_dataset,
    resolve_dataset_id,
    species_configs,
    species_id_for_dataset,
    variants_for_species,
)
from keyword_search import (
    cross_species_best_annotations_for_genes,
    keyword_index_available,
    search_keywords,
)

from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from flask import (
    Flask,
    Response,
    jsonify,
    redirect,
    render_template,
    request,
    session,
    url_for,
)

app = Flask(__name__)
app.secret_key = os.environ.get("FLASK_SECRET_KEY", "panviewer-dev-session-key")

# Latest FAMSA MSA per browser session: session holds a token; full alignment lives
# here (process-local). If the worker has no entry (e.g. another Gunicorn worker
# served the page), ``/api/newick_fasttree`` re-runs FAMSA for the requested gene IDs.
_FAMSA_MSA_MEM: OrderedDict[str, tuple[dict[str, str], float]] = OrderedDict()
_FAMSA_MSA_MEM_LOCK = threading.Lock()
_FAMSA_MSA_MEM_MAX = 80
_FAMSA_MSA_MEM_TTL_SEC = 24 * 3600.0


def _evict_famsa_msa_mem_unlocked(now: float) -> None:
    dead = [
        cid
        for cid, (_, t0) in _FAMSA_MSA_MEM.items()
        if now - t0 > _FAMSA_MSA_MEM_TTL_SEC
    ]
    for cid in dead:
        del _FAMSA_MSA_MEM[cid]
    while len(_FAMSA_MSA_MEM) > _FAMSA_MSA_MEM_MAX:
        _FAMSA_MSA_MEM.popitem(last=False)


def register_latest_famsa_alignment(aligned: dict[str, str]) -> None:
    """Store the latest full MSA in RAM and record its id in the Flask session."""
    if len(aligned) < 2:
        session.pop("famsa_msa_id", None)
        session.modified = True
        return
    snap = {str(k): str(v) for k, v in aligned.items()}
    cid = secrets.token_urlsafe(24)
    now = time.time()
    with _FAMSA_MSA_MEM_LOCK:
        _evict_famsa_msa_mem_unlocked(now)
        while len(_FAMSA_MSA_MEM) >= _FAMSA_MSA_MEM_MAX:
            _FAMSA_MSA_MEM.popitem(last=False)
        _FAMSA_MSA_MEM[cid] = (snap, now)
        _FAMSA_MSA_MEM.move_to_end(cid)
    session["famsa_msa_id"] = cid
    session.modified = True


def get_latest_famsa_alignment_from_session() -> dict[str, str] | None:
    """Return a copy of the MSA last registered for this session, if still cached."""
    cid = session.get("famsa_msa_id")
    if not isinstance(cid, str):
        return None
    now = time.time()
    with _FAMSA_MSA_MEM_LOCK:
        _evict_famsa_msa_mem_unlocked(now)
        hit = _FAMSA_MSA_MEM.get(cid)
        if not hit:
            return None
        snap, t0 = hit
        if now - t0 > _FAMSA_MSA_MEM_TTL_SEC:
            try:
                del _FAMSA_MSA_MEM[cid]
            except KeyError:
                pass
            return None
        return dict(snap)


def famsa_rebuild_and_cache_for_genes(
    dataset_id: str, gene_ids: list[str]
) -> dict[str, str] | None:
    """Re-run FAMSA for these gene IDs, register the MSA in session RAM, return full aligned dict."""
    seen: set[str] = set()
    clean: list[str] = []
    for g in gene_ids:
        s = str(g).strip()
        if s and s not in seen:
            seen.add(s)
            clean.append(s)
    if len(clean) < 2 or len(clean) > 600:
        return None
    conn = get_db(dataset_id)
    try:
        merged, _acc = collect_sequences_for_gene_ids(conn, clean, dataset_id)
    finally:
        conn.close()
    if len(merged) < 2:
        return None
    aligned, _nw = run_famsa(merged)
    if len(aligned) < 2:
        return None
    register_latest_famsa_alignment(aligned)
    return aligned


def wheat_pan_display_label(pan_id: str | None) -> str:
    """Wheat UI label; DB ids are already ``Traes_pan*`` (legacy ``pan_*`` URLs still accepted)."""
    s = (pan_id or "").strip()
    low = s.lower()
    if low.startswith("traes_"):
        return s
    if len(s) >= 4 and low.startswith("pan_"):
        return "Traes_pan" + s[4:]
    if low.startswith("pan"):
        return f"Traes_{s}"
    return s


def oat_pan_display_label(pan_id: str | None) -> str:
    """Oat Pandagma UI label; DB ids are ``Avena_pan*`` (PanOat uses ``Avena_N0.HOG*``)."""
    s = (pan_id or "").strip()
    low = s.lower()
    if low.startswith("avena_"):
        return s
    if len(s) >= 4 and low.startswith("pan_"):
        return "Avena_pan" + s[4:]
    return s


def internal_pangene_id(pangene_id: str | None, dataset_id: str | None) -> str:
    """Map route pan IDs to the form stored in ``genes.pangene`` for this variant."""
    s = (pangene_id or "").strip()
    if not s:
        return s
    sid = species_id_for_dataset((dataset_id or "").strip().lower())
    low = s.lower()
    if sid == "wheat":
        if low.startswith("traes_"):
            return s
        if low.startswith("pan"):
            return wheat_pan_display_label(s)
    if sid == "oat":
        if low.startswith("avena_"):
            return s
        if low.startswith("pan"):
            return oat_pan_display_label(s)
    if sid == "barley":
        if low.startswith("horvu_"):
            return s
        if low.startswith("pan"):
            pl = s.lower()
            if pl.startswith("pan_"):
                return "HORVU_pan" + s[4:]
            return f"HORVU_{s}"
    return s


def pangene_display_label(pan_id: str | None, dataset_id: str | None) -> str:
    """Species-specific pan-gene label for titles, breadcrumbs, and search suggestions."""
    sid = species_id_for_dataset((dataset_id or "").strip().lower())
    if sid == "wheat":
        return wheat_pan_display_label(pan_id)
    if sid == "oat":
        return oat_pan_display_label(pan_id)
    return (pan_id or "").strip()


app.jinja_env.globals["wheat_pan_display"] = wheat_pan_display_label
app.jinja_env.globals["oat_pan_display"] = oat_pan_display_label
app.jinja_env.globals["pangene_display"] = pangene_display_label


def app_root_url_for(endpoint: str, **values: Any) -> str:
    """
    Like :func:`flask.url_for`, but if the path omits :data:`APPLICATION_ROOT` (some proxies
    clear ``SCRIPT_NAME``), prepend it. Otherwise ``fetch`` posts to ``/api/...`` at the
    site root and can hit the wrong vhost (405 Method Not Allowed on POST).
    """
    p = str(url_for(endpoint, **values))
    root = (APPLICATION_ROOT or "").rstrip("/")
    if not root or p.startswith(root + "/") or p == root:
        return p
    if p.startswith("/") and not p.startswith("//"):
        return f"{root}{p}"
    return p


app.jinja_env.globals["app_root_url_for"] = app_root_url_for


def grain_genes_blast_transfer_api_url() -> str:
    """Path to PanViewer proxy for GrainGenes ``query_transfer`` (no url_for)."""
    root = (APPLICATION_ROOT or "").rstrip("/")
    return f"{root}/api/grain_genes_blast_transfer"


app.jinja_env.globals["grain_genes_blast_transfer_api_url"] = grain_genes_blast_transfer_api_url


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
DATABASE_DIR = os.path.join(BASE_DIR, "database")
INPUT_DIR = os.path.join(BASE_DIR, "input")
DEFAULT_DATASET_ID = "wheat"


def genome_subgenome_layout(dataset_id: str) -> str:
    """Gene grid: wheat A/B/D; oat hexaploid A/C/D; other species use a single Chr column."""
    sid = species_id_for_dataset((dataset_id or "").strip().lower())
    if sid == "wheat":
        return "abd"
    if sid == "oat":
        return "acd"
    return "chr"


def _subgenome_letter_from_chr_name(chrom: str | None) -> str | None:
    """
    Parse trailing or leading A/B/C/D from chromosome names (e.g. ``3C``, ``chr1A``, ``D4``).
    Used when ``gene_coords.subgenome`` is unset but ``chr`` encodes the subgenome arm.
    """
    if chrom is None:
        return None
    s = str(chrom).strip().upper()
    if not s:
        return None
    if s.startswith("CHR"):
        s = s[3:]
    m = re.fullmatch(r"\d+([ABCD])", s)
    if m:
        return m.group(1)
    m = re.fullmatch(r"([ABCD])\d+", s)
    if m:
        return m.group(1)
    return None


def subgenome_for_ui_from_coords_row(cr: dict | None, dataset_id: str) -> str | None:
    """A/B/D or A/C/D letter for grids, using DB ``subgenome`` or inferring from ``chr``."""
    if not cr:
        return None
    sg = (cr.get("subgenome") or "").strip().upper()
    if not sg:
        inferred = _subgenome_letter_from_chr_name(cr.get("chr"))
        if inferred:
            sg = inferred
    layout = genome_subgenome_layout(dataset_id)
    if layout == "abd":
        return sg if sg in ("A", "B", "D") else None
    if layout == "acd":
        return sg if sg in ("A", "C", "D") else None
    return None


def load_external_db_links() -> list[dict[str, str]]:
    """
    Rows from static/links.csv for external DB buttons on the Genes tab.
    URLs may contain {gene_id} and/or {chr}, {start}, {end} placeholders.
    Optional ``genome`` column (wheat, barley, oat) scopes rows to that dataset;
    rows with an empty genome match all datasets (legacy CSV).
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
                        "genome": (row.get("genome") or "").strip().lower(),
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


def external_db_links_for_dataset(dataset_id: str) -> list[dict[str, str]]:
    """Subset of ``EXTERNAL_DB_LINKS`` for the current species (see ``genome`` in links.csv)."""
    ds = species_id_for_dataset((dataset_id or "").strip().lower())
    rows: list[dict[str, str]] = []
    for row in EXTERNAL_DB_LINKS:
        genome = (row.get("genome") or "").strip().lower()
        if not genome or genome == ds:
            rows.append(row)
    return rows


def load_genome_metadata() -> list[dict[str, str]]:
    """
    Rows from static/genomes.csv for the Accessions tab (source, version, link).
    Optional ``genome`` column (wheat, barley, oat) scopes rows to that species.
    """
    path = os.path.join(BASE_DIR, "static", "genomes.csv")
    out: list[dict[str, str]] = []
    if not os.path.isfile(path):
        return out
    try:
        with open(path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for row in reader:
                accession = (row.get("accession") or "").strip()
                if not accession:
                    continue
                out.append(
                    {
                        "genome": (row.get("genome") or "").strip().lower(),
                        "accession": accession,
                        "source": (row.get("source") or "").strip(),
                        "version": (row.get("version") or "").strip(),
                        "link": (row.get("link") or "").strip(),
                    }
                )
    except OSError:
        return out
    return out


GENOME_METADATA = load_genome_metadata()


def genome_metadata_for_dataset(dataset_id: str) -> list[dict[str, str]]:
    """Subset of ``GENOME_METADATA`` for the current species (see ``genome`` in genomes.csv)."""
    ds = species_id_for_dataset((dataset_id or "").strip().lower())
    rows: list[dict[str, str]] = []
    for row in GENOME_METADATA:
        genome = (row.get("genome") or "").strip().lower()
        if not genome or genome == ds:
            rows.append(row)
    return rows


def grain_genes_blast_config_for_dataset(dataset_id: str) -> dict[str, str]:
    """
    GrainGenes BLAST ``query_transfer`` database ids for the active dataset.

    IDs match SequenceServer categories on https://graingenes.org/blast/ (see
    gg-db-groups.json / searchdata.json titles).
    """
    did = (dataset_id or "").strip().lower()
    base = "https://graingenes.org/blast/"
    if did == "panoat":
        return {
            "base_url": base,
            "nucl_db": "AvSangChr",
            "prot_db": "AvSangChr",
            "nucl_label": "PanOat Sang v1.1",
            "prot_label": "PanOat Sang v1.1",
        }
    sid = species_id_for_dataset(did)
    if sid == "wheat":
        return {
            "base_url": base,
            "nucl_db": "IWGSCv2",
            "prot_db": "IWGSCv2Prot",
            "nucl_label": "Chinese Spring IWGSC RefSeq v2.1",
            "prot_label": "Chinese Spring IWGSC RefSeq v2.1 proteins",
        }
    if sid == "barley":
        return {
            "base_url": base,
            "nucl_db": "Hv-Morex3",
            "prot_db": "Hv-Morex3-Prot",
            "nucl_label": "Morex v3",
            "prot_label": "Morex v3 proteins",
        }
    if sid == "oat":
        return {
            "base_url": base,
            "nucl_db": "Asativa-sang",
            "prot_db": "Asativa-sang",
            "nucl_label": "Sang v1.1",
            "prot_label": "Sang v1.1",
        }
    return {"base_url": base, "nucl_db": "", "prot_db": "", "nucl_label": "", "prot_label": ""}


GRAINGENES_BLAST_TRANSFER_URL = "https://graingenes.org/blast/query_transfer"


def load_about_stats(dataset_id: str) -> dict[str, str] | None:
    """
    Read ``database/stats.tsv`` (one row per species from ``build_index.py``).

    Expects columns ``accessions``, ``pan_genes``, ``genes``.
    Single-row files without ``species`` are treated as wheat.
    """
    dataset_id = (dataset_id or "").strip().lower()
    path = os.path.join(DATABASE_DIR, "stats.tsv")
    if not os.path.isfile(path):
        return None

    def norm_row(r: dict) -> dict[str, str]:
        return {((k or "").strip().lower()): (v or "") for k, v in r.items()}

    def coalesce_stats(nr: dict[str, str]) -> dict[str, str] | None:
        acc = str(nr.get("accessions") or "").strip()
        genes = str(nr.get("genes") or "").strip()
        pan = str(nr.get("pan_genes") or "").strip()
        if not all((acc, pan, genes)):
            return None
        try:
            return {
                "accessions": f"{int(acc):,}",
                "pan_genes": f"{int(pan):,}",
                "genes": f"{int(genes):,}",
            }
        except ValueError:
            return None

    try:
        with open(path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            fieldnames = [x.strip().lower() for x in (reader.fieldnames or [])]
            if "species" in fieldnames:
                for row in reader:
                    nr = norm_row(row)
                    sid = (nr.get("species") or "").strip().lower()
                    if sid != dataset_id:
                        continue
                    return coalesce_stats(nr)
                return None
            if dataset_id != "wheat":
                return None
            row = next(reader, None)
            if not row:
                return None
            return coalesce_stats(norm_row(row))
    except (OSError, ValueError, TypeError):
        return None


def available_datasets():
    out = set()
    for ds_id, cfg in dataset_configs().items():
        if not os.path.exists(cfg.get("db_path", "")):
            continue

        mode = cfg.get("mode")
        if mode == "pandagma":
            # Sequences and synteny are served from the DB (e.g. protein_seq_map, gene_coords).
            # FASTA under input/<species>/prot is only for index builds, not runtime.
            out.add(ds_id)
            continue

        # N0 + FASTA directory mode requires a per-cluster protein FASTA directory.
        fasta_dir = cfg.get("fasta_dir")
        if fasta_dir and os.path.isdir(fasta_dir):
            out.add(ds_id)
    return out


def get_dataset_config(dataset_id: str):
    dataset_id = resolve_dataset_id((dataset_id or "").strip().lower()) or ""
    cfg = dataset_configs()
    if dataset_id in cfg:
        return cfg[dataset_id]
    if DEFAULT_DATASET_ID in cfg:
        return cfg[DEFAULT_DATASET_ID]
    if cfg:
        first = sorted(cfg.keys())[0]
        return cfg[first]
    raise RuntimeError("No SQLite databases found under database/")


def _open_sqlite_connection(path: str) -> sqlite3.Connection:
    """
    Open a species database read-only (runtime never writes).

    Uses URI ``mode=ro&immutable=1&nolock=1`` so SQLite does not create or update
    ``-wal`` / ``-shm`` sidecars — safe for Docker ``:ro`` bind mounts.
    """
    path = os.path.abspath(path)
    uri = f"file:{path}?mode=ro&immutable=1&nolock=1"
    last_err: sqlite3.OperationalError | None = None
    for attempt in range(3):
        try:
            conn = sqlite3.connect(uri, uri=True, timeout=30.0)
            conn.row_factory = sqlite3.Row
            conn.execute("PRAGMA busy_timeout=30000")
            conn.execute("SELECT 1")
            return conn
        except sqlite3.OperationalError as e:
            last_err = e
            time.sleep(0.05 * (attempt + 1))
    raise RuntimeError(f"Cannot open database {path}: {last_err}") from last_err


def get_db(dataset_id: str):
    cfg = get_dataset_config(dataset_id)
    path = os.path.abspath(cfg["db_path"])
    if not os.path.isfile(path):
        raise RuntimeError(
            f"Database file not found for dataset '{dataset_id}': {path}"
        )
    return _open_sqlite_connection(path)


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


_genes_pangene_ok_lock = threading.Lock()
_genes_pangene_ok_by_db_path: set[str] = set()

_pangene_info_meta_lock = threading.Lock()
_pangene_info_meta_by_db_path: dict[str, tuple[str, str] | None] = {}


def ensure_genes_has_pangene(conn: sqlite3.Connection) -> None:
    """Require ``genes.pangene`` (and optional ``pangene_info``) from ``build_index.py``."""
    key = _main_sqlite_db_path(conn)
    with _genes_pangene_ok_lock:
        if key in _genes_pangene_ok_by_db_path:
            return
        cols = {r[1] for r in conn.execute("PRAGMA table_info(genes)").fetchall()}
        if "pan_gene" in cols and "pangene" not in cols:
            raise RuntimeError(
                f"SQLite genes table at {key} uses legacy column 'pan_gene'. "
                "Rebuild with build_index.py."
            )
        if "pangene" not in cols:
            raise RuntimeError(
                f"SQLite genes table at {key} must include column 'pangene'. "
                "Rebuild with build_index.py."
            )
        bad = conn.execute(
            "SELECT 1 FROM genes WHERE pangene IS NULL OR TRIM(pangene) = '' LIMIT 1"
        ).fetchone()
        if bad:
            raise RuntimeError(
                f"SQLite genes table at {key} has rows with empty pangene. "
                "Rebuild with build_index.py."
            )
        has_pi = conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='pangene_info' LIMIT 1"
        ).fetchone()
        if has_pi:
            pic_pi = {r[1] for r in conn.execute("PRAGMA table_info(pangene_info)").fetchall()}
            if "pangene" not in pic_pi:
                raise RuntimeError(
                    f"SQLite pangene_info at {key} exists but has no 'pangene' column. "
                    "Rebuild with build_index.py."
                )
            bad_pi = conn.execute(
                "SELECT 1 FROM pangene_info WHERE pangene IS NULL OR TRIM(pangene) = '' LIMIT 1"
            ).fetchone()
            if bad_pi:
                raise RuntimeError(
                    f"SQLite pangene_info at {key} has rows with empty pangene. "
                    "Rebuild with build_index.py."
                )
        _genes_pangene_ok_by_db_path.add(key)


def pangene_info_meta(conn: sqlite3.Connection) -> tuple[str, str] | None:
    """If ``pangene_info`` exists with a ``pangene`` column, return (table, id_column)."""
    key = _main_sqlite_db_path(conn)
    with _pangene_info_meta_lock:
        hit = _pangene_info_meta_by_db_path.get(key)
    if hit is not None:
        return hit
    meta: tuple[str, str] | None = None
    row = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='pangene_info' LIMIT 1"
    ).fetchone()
    if row:
        pic = {r[1] for r in conn.execute("PRAGMA table_info(pangene_info)").fetchall()}
        if "pangene" in pic:
            meta = ("pangene_info", "pangene")
    with _pangene_info_meta_lock:
        _pangene_info_meta_by_db_path[key] = meta
    return meta


def extract_accession(gene_id):
    """Extract accession name from gene ID like HORVU.BONUS.PROJ.1HG00046750.1"""
    parts = gene_id.split(".")
    if len(parts) >= 3:
        return parts[1]
    return gene_id


def escape_sql_like(pattern: str) -> str:
    """Escape ``%`` and ``_`` for SQLite ``LIKE ... ESCAPE '\\'`` (underscore matches one char by default)."""
    return (
        (pattern or "")
        .replace("\\", "\\\\")
        .replace("%", "\\%")
        .replace("_", "\\_")
    )


def gene_id_display_label(gene_id: str) -> str:
    """Prefer the transcript-style id after ``accession|`` for compact UI (barley / pandagma BED)."""
    s = (gene_id or "").strip()
    if "|" in s:
        return s.split("|", 1)[1].strip() or s
    return s


def morex_v3_reference_gene_ids_from_rows(genes) -> list[str]:
    """Full DB ``gene_id`` values for MorexV3 (reference) rows — used for Expression + alignment keys."""
    seen: set[str] = set()
    out: list[str] = []
    for row in genes:
        acc = (row["accession"] or "").strip()
        if acc.lower() != "morexv3":
            continue
        gid = (row["gene_id"] or "").strip()
        if gid and gid not in seen:
            seen.add(gid)
            out.append(gid)
    out.sort()
    return out


def chinese_spring_gene_ids_from_rows(genes) -> list[str]:
    """All genes in this pan-gene cluster with Chinese Spring accession (no space), case-insensitive, sorted."""
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


def sang_v11_reference_gene_ids_from_rows(genes) -> list[str]:
    """Gene IDs for sangV11 reference rows (oat expression via PlantApp genome AsSang)."""
    seen: set[str] = set()
    out: list[str] = []
    for row in genes:
        acc = re.sub(r"[^a-z0-9]", "", (row["accession"] or "").strip().lower())
        if acc != "sangv11":
            continue
        gid = (row["gene_id"] or "").strip()
        if gid and gid not in seen:
            seen.add(gid)
            out.append(gid)
    out.sort()
    return out


def _looks_like_barley_gene_query(q: str) -> bool:
    """Heuristic: barley Morex-style locus ids start with ``HORVU.``."""
    s = (q or "").strip()
    return bool(re.match(r"^HORVU\.", s, re.I))


def _looks_like_wheat_gene_query(q: str) -> bool:
    """Heuristic: wheat IWGSC-style ids start with ``TraesCS`` (any case)."""
    s = (q or "").strip()
    return bool(re.match(r"^TraesCS", s, re.I))


def _looks_like_oat_gene_query(q: str) -> bool:
    """Heuristic: oat AVESA-style locus ids start with ``AVESA.``."""
    s = (q or "").strip()
    return bool(re.match(r"^AVESA\.", s, re.I))


_UNSAFE_QUERY_CHARS_RE = re.compile(r"""[;:{}\(\)\*<>'"`\\|?%#=\[\]^$@!~/&+]""")


def strip_unsafe_query_chars(query: str) -> str:
    return _UNSAFE_QUERY_CHARS_RE.sub("", query or "").strip()


def search_genes(query, dataset_id: str):
    """Search for genes matching the query (exact or partial)."""
    q = (query or "").strip()
    conn = get_db(dataset_id)
    ensure_genes_has_pangene(conn)
    cur = conn.cursor()
    cur.execute(
        f"SELECT gene_id, pangene AS pan, accession FROM genes WHERE gene_id = ? COLLATE NOCASE",
        (canonical_gene_id(q),),
    )
    results = cur.fetchall()

    if not results and q:
        like_pat = f"%{escape_sql_like(canonical_gene_id(q))}%"
        cur.execute(
            f"SELECT gene_id, pangene AS pan, accession FROM genes "
            "WHERE gene_id LIKE ? ESCAPE '\\' COLLATE NOCASE LIMIT 200",
            (like_pat,),
        )
        results = cur.fetchall()

    conn.close()
    return results


def get_genes_for_pangene(pangene_id: str, dataset_id: str):
    """
    Member genes for one pan-gene cluster.

    Pandagma SQLite stores the pan-gene id in ``genes.pangene``.
    """
    conn = get_db(dataset_id)
    ensure_genes_has_pangene(conn)
    cur = conn.cursor()
    if table_exists(conn, "gene_coords"):
        cur.execute(
            f"""
            SELECT g.gene_id AS gene_id, g.accession AS accession,
                   gc.chr AS chr, gc.start AS start, gc.end AS end
            FROM genes g
            LEFT JOIN gene_coords gc
              ON g.gene_id = gc.gene_id AND g.accession = gc.accession
            WHERE g.pangene = ?
            ORDER BY g.accession
            """,
            (pangene_id,),
        )
    else:
        cur.execute(
            f"SELECT gene_id, accession FROM genes WHERE pangene = ? ORDER BY accession",
            (pangene_id,),
        )
    results = cur.fetchall()
    conn.close()
    return results


def build_wheat_cluster_picker_payload(
    conn,
    scope_pangene: str,
    current_pangene_id: str,
    *,
    max_pangenes: int = 200,
) -> dict | None:
    """
    Data for wheat pan-gene detail: genes in ``scope_pangene`` (one pan-gene id), with optional A/B/D dominants.
    Used for multi-select UI (accessions/genes/synteny scope). ``scope_pangene`` is the page pan-gene id.
    """
    ensure_genes_has_pangene(conn)
    cur = conn.cursor()
    cur.execute(
        "SELECT COUNT(DISTINCT pangene) AS n FROM genes WHERE pangene = ?", (scope_pangene,)
    )
    nh = int(cur.fetchone()["n"])
    if nh == 0:
        return None
    if nh > max_pangenes:
        return {
            "pangene": scope_pangene,
            "skipped": True,
            "reason": (
                f"This pan-gene scope has {nh} distinct pan-gene rows; the picker supports up to {max_pangenes}."
            ),
            "pangene_count": nh,
            "current_pangene": current_pangene_id,
        }

    coords_on = table_exists(conn, "gene_coords")
    cur.execute(
        "SELECT pangene AS pan, gene_id, accession FROM genes WHERE pangene = ? ORDER BY pangene, gene_id",
        (scope_pangene,),
    )
    by_pan: dict[str, list[dict]] = defaultdict(list)
    for r in cur.fetchall():
        pan_id = (r["pan"] or "").strip()
        ge: dict = {
            "gene_id": r["gene_id"],
            "accession": r["accession"],
            "pangene": pan_id,
        }
        if coords_on:
            summ = _gene_coords_summary(cur, ge["gene_id"], ge["accession"])
            ge["coords"] = summ
            cr = fetch_coords_row(cur, ge["gene_id"], ge["accession"])
            if cr:
                ge["subgenome"] = subgenome_for_ui_from_coords_row(cr, "wheat")
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
        by_pan[pan_id].append(ge)
    pans_sorted = sorted(by_pan.keys())

    panel = None
    if coords_on:
        panel = build_wheat_cluster_pangene_triad_panel(
            conn, scope_pangene, current_pangene_id, max_pangenes=max_pangenes
        )
    row_by_pangene: dict[str, dict] = {}
    if panel and not panel.get("skipped"):
        for r in panel["rows"]:
            row_by_pangene[r["pangene"]] = r

    items: list[dict] = []
    for h in pans_sorted:
        pr = row_by_pangene.get(h)
        gl = by_pan[h]
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
        "pangene": scope_pangene,
        "skipped": False,
        "pangene_items": items,
        "current_pangene": current_pangene_id,
    }


LARGE_PANGENE_GENE_THRESHOLD = 500


def count_pangene_genes(pangene_id: str, dataset_id: str) -> int:
    """Member count for one pan-gene (``pangene_info.gene_count`` when indexed)."""
    conn = get_db(dataset_id)
    try:
        ensure_genes_has_pangene(conn)
        cur = conn.cursor()
        meta = pangene_info_meta(conn)
        if meta:
            tbl, col = meta
            cur.execute(f"SELECT gene_count FROM {tbl} WHERE {col} = ?", (pangene_id,))
            row = cur.fetchone()
            if row and row[0] is not None:
                return int(row[0])
        cur.execute(
            "SELECT COUNT(*) FROM genes WHERE pangene = ?",
            (pangene_id,),
        )
        return int(cur.fetchone()[0])
    finally:
        conn.close()


def get_pangene_info(pangene_id: str, dataset_id: str):
    """
    Return metadata for one pan-gene cluster.

    Uses ``pangene_info`` when present; otherwise checks membership in ``genes``.
    """
    conn = get_db(dataset_id)
    try:
        ensure_genes_has_pangene(conn)
        cur = conn.cursor()
        meta = pangene_info_meta(conn)
        if meta:
            tbl, col = meta
            cur.execute(f"SELECT * FROM {tbl} WHERE {col} = ?", (pangene_id,))
            r = cur.fetchone()
            if r:
                return r
        cur.execute(
            "SELECT 1 FROM genes WHERE pangene = ? COLLATE NOCASE LIMIT 1",
            (pangene_id,),
        )
        if cur.fetchone():
            return {"pangene": pangene_id}
        return None
    finally:
        conn.close()


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


def read_fasta(pangene_id, dataset_id: str):
    """
    Read protein sequences for a given group id (Pandagma pan_id in Pandagma mode).
    """
    cfg = get_dataset_config(dataset_id)

    fasta_dir = cfg.get("fasta_dir")
    fasta_path = (
        os.path.join(fasta_dir, f"{pangene_id}.protein.fasta")
        if fasta_dir
        else None
    )
    if fasta_path and os.path.exists(fasta_path):
        sequences = {}
        for record in SeqIO.parse(fasta_path, "fasta"):
            gid = None
            for tok in fasta_header_token_candidates(record):
                c = canonical_gene_id(tok)
                if c:
                    gid = c
                    break
            if not gid:
                continue
            sequences[gid] = str(record.seq)
        return sequences

    if cfg.get("mode") == "pandagma":
        conn = get_db(dataset_id)
        try:
            if not table_exists(conn, "protein_seq_map"):
                return {}
            ensure_genes_has_pangene(conn)
            cur = conn.cursor()
            cur.execute(
                f"""
                SELECT g.gene_id, u.seq_comp
                FROM genes g
                JOIN protein_seq_map m ON m.gene_id = g.gene_id
                JOIN protein_seq_uniq u ON u.uniq_id = m.uniq_id
                WHERE g.pangene = ?
                ORDER BY g.accession, g.gene_id
                """,
                (pangene_id,),
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
    Load ungapped protein sequences for the given gene IDs (any pan-gene cluster in the DB).
    Returns (gene_id -> sequence, gene_id -> accession).
    """
    cur = conn.cursor()
    ensure_genes_has_pangene(conn)
    pan_to_ids: dict[str, set[str]] = {}
    acc_map: dict[str, str] = {}
    for gid in gene_ids:
        row = None
        for vid in gene_id_coord_variants(gid):
            cur.execute(
                f"""
                SELECT gene_id, pangene AS pan, accession FROM genes
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
        pan = (row["pan"] or "").strip()
        acc = (row["accession"] or "").strip()
        if not db_g or not pan:
            continue
        pan_to_ids.setdefault(pan, set()).add(db_g)
        acc_map[db_g] = acc or db_g

    merged: dict[str, str] = {}
    for pan, idset in pan_to_ids.items():
        seqs = read_fasta(pan, dataset_id)
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
# FastTree on FAMSA MSA (seconds).
FASTTREE_TIMEOUT_SEC = 180


def _find_fasttree_bin() -> str | None:
    return shutil.which("FastTreeMP") or shutil.which("FastTree") or shutil.which("fasttree")


def _newick_output_file_valid(path: str) -> bool:
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return False
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            s = f.read().lstrip()
    except OSError:
        return False
    return bool(s) and s[0] == "("


def _fasttree_newick_from_fasta(
    fasta_path: str, newick_path: str, log: logging.Logger, timeout: int
) -> bool:
    """Run FastTree ``-lg -quiet`` on a protein FASTA (fallback: FAMSA MSA)."""
    ft = _find_fasttree_bin()
    if not ft:
        log.warning("FastTree not on PATH; conda install -c bioconda fasttree")
        return False
    if os.path.isfile(newick_path):
        try:
            os.unlink(newick_path)
        except OSError:
            pass
    # Protein MSA: -lg (Le+Gascuel); -quiet avoids stderr chatter.
    cmd = [ft, "-lg", "-quiet", fasta_path]
    try:
        with open(newick_path, "w", encoding="utf-8") as out:
            r = subprocess.run(
                cmd,
                stdout=out,
                stderr=subprocess.PIPE,
                text=True,
                timeout=timeout,
            )
    except OSError as exc:
        log.warning("FastTree subprocess error: %s", exc)
        return False
    if r.returncode != 0 or not _newick_output_file_valid(newick_path):
        err = (r.stderr or r.stdout or "").strip()[:800]
        log.warning(
            "FastTree failed on %s (exit %s): %s",
            os.path.basename(fasta_path),
            r.returncode,
            err or "(no stderr)",
        )
        return False
    return True


def run_famsa(sequences, *, timeout_sec=FAMSA_TIMEOUT_SEC):
    """Align with FAMSA using neighbour-joining guide export/import when supported.

    Primary path (FAMSA2): ``-gt nj -gt_export`` to a temp Newick, then
    ``-gt import`` into the alignment. The returned tree is that NJ Newick after
    :func:`midpoint_root_newick`. If export/import is unavailable or fails,
    falls back to default FAMSA alignment plus **FastTree** on the MSA (also
    midpoint-rooted).
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

    tmp_nj = tmp_in_path + ".famsa_nj.nwk"
    tmp_aln = tmp_in_path + ".famsa.aln"
    tmp_aln_tree = tmp_in_path + ".from_aln.nwk"

    def _cleanup():
        for p in (tmp_in_path, tmp_nj, tmp_aln, tmp_aln_tree):
            if os.path.exists(p):
                try:
                    os.unlink(p)
                except OSError:
                    pass

    # Do not pass -keep_duplicates: FAMSA 1.x ignores unknown flags and treats the
    # next token as the input path, yielding "Unable to open input file -keep_duplicates".
    famsa_export_nj = [
        famsa_bin,
        "-t",
        "0",
        "-gt",
        "nj",
        "-gt_export",
        tmp_in_path,
        tmp_nj,
    ]
    famsa_import = [
        famsa_bin,
        "-t",
        "0",
        "-gt",
        "import",
        tmp_nj,
        tmp_in_path,
        tmp_aln,
    ]
    famsa_default = [famsa_bin, "-t", "0", tmp_in_path, tmp_aln]

    def _run_famsa_nj_then_midpoint() -> tuple[dict[str, str], str] | None:
        try:
            r0 = subprocess.run(
                famsa_export_nj,
                capture_output=True,
                text=True,
                timeout=timeout_sec,
            )
        except FileNotFoundError:
            return None
        if r0.returncode != 0 or not _newick_output_file_valid(tmp_nj):
            tail = (r0.stderr or r0.stdout or "").strip()[:400]
            log.debug(
                "famsa NJ export not used (%s); falling back to default + FastTree",
                tail or f"exit {r0.returncode}",
            )
            return None
        try:
            r2 = subprocess.run(
                famsa_import,
                capture_output=True,
                text=True,
                timeout=timeout_sec,
            )
        except FileNotFoundError:
            return None
        if r2.returncode != 0 or not os.path.isfile(tmp_aln):
            err = (r2.stderr or r2.stdout or "").strip()[:400]
            log.debug("famsa -gt import failed (%s); fallback", err or f"exit {r2.returncode}")
            return None
        aligned_local: dict[str, str] = {}
        for record in SeqIO.parse(tmp_aln, "fasta"):
            aligned_local[record.id] = str(record.seq)
        if not aligned_local:
            return None
        try:
            with open(tmp_nj, encoding="utf-8", errors="replace") as f:
                raw_nj = f.read().strip()
        except OSError:
            raw_nj = ""
        if not raw_nj:
            return None
        return aligned_local, midpoint_root_newick(raw_nj)

    def _run_default_plus_fasttree() -> tuple[dict[str, str], str]:
        newick_ft = ""
        try:
            r2 = subprocess.run(
                famsa_default,
                capture_output=True,
                text=True,
                timeout=timeout_sec,
            )
        except FileNotFoundError:
            return sequences, ""
        if r2.returncode != 0 or not os.path.isfile(tmp_aln):
            err = (r2.stderr or r2.stdout or "").strip()[:800]
            log.warning("famsa alignment failed: %s", err or "(no output)")
            return sequences, ""
        aligned_local: dict[str, str] = {}
        for record in SeqIO.parse(tmp_aln, "fasta"):
            aligned_local[record.id] = str(record.seq)
        if not aligned_local:
            return sequences, ""
        if _fasttree_newick_from_fasta(tmp_aln, tmp_aln_tree, log, FASTTREE_TIMEOUT_SEC):
            with open(tmp_aln_tree, encoding="utf-8", errors="replace") as f:
                newick_ft = f.read().strip()
        else:
            log.warning(
                "FastTree on FAMSA MSA failed; alignment returned without a tree"
            )
        if newick_ft:
            newick_ft = midpoint_root_newick(newick_ft)
        return aligned_local, newick_ft

    try:
        hit = _run_famsa_nj_then_midpoint()
        if hit is not None:
            return hit
        return _run_default_plus_fasttree()
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


def midpoint_root_newick(newick_str: str) -> str:
    """Midpoint-root a Newick tree (no outgroup); no-op on parse failure.

    Used for FAMSA NJ and FastTree outputs so the UI root is arbitrary but
    balanced, not an outlier taxon.
    """
    s = (newick_str or "").strip()
    if not s:
        return ""
    try:
        tree = Phylo.read(io.StringIO(s), "newick")
        if len(tree.get_terminals()) < 2:
            return s
        rm = getattr(tree, "root_at_midpoint", None)
        if rm is None:
            return s
        rm()
        buf = io.StringIO()
        Phylo.write(tree, buf, "newick")
        out = buf.getvalue().strip()
        line = out.splitlines()[0].strip() if out else ""
        return line if line else s
    except Exception:
        return s


def get_tree_leaf_order(newick_str):
    """Parse a Newick tree and return leaf names in display order (top to bottom)."""
    if not newick_str:
        return []
    try:
        tree = Phylo.read(io.StringIO(newick_str), "newick")
        return [clade.name for clade in tree.get_terminals()]
    except Exception:
        return []


def _sanitize_aligned_for_fasttree_api(raw: dict) -> dict[str, str]:
    """Keep only plausible gene_id → aligned protein strings from JSON."""
    out: dict[str, str] = {}
    if not isinstance(raw, dict):
        return out
    for k, v in raw.items():
        gid = str(k).strip()
        if not gid or len(gid) > 256:
            continue
        s = str(v).strip()
        if not s or len(s) > 120_000:
            continue
        if not re.fullmatch(r"[A-Za-z\-.*]+", s):
            continue
        out[gid] = s
    return out


def newick_fasttree_from_aligned(
    aligned: dict[str, str],
    *,
    timeout: int = FASTTREE_TIMEOUT_SEC,
) -> tuple[str, list[str]]:
    """Run FastTree on an existing MSA; midpoint-root; return (newick, leaf_order)."""
    log = logging.getLogger(__name__)
    if len(aligned) < 2:
        return "", []
    tmp_fa = ""
    tmp_nwk = ""
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".msa.fa", delete=False, encoding="utf-8"
        ) as tmp:
            for gid in sorted(aligned.keys()):
                tmp.write(f">{gid}\n{aligned[gid]}\n")
            tmp_fa = tmp.name
        tmp_nwk = tmp_fa + ".fasttree.nwk"
        if not _fasttree_newick_from_fasta(tmp_fa, tmp_nwk, log, timeout):
            return "", []
        with open(tmp_nwk, encoding="utf-8", errors="replace") as f:
            raw = f.read().strip()
        if not raw:
            return "", []
        nw = midpoint_root_newick(raw)
        return nw, get_tree_leaf_order(nw)
    except OSError as exc:
        log.warning("newick_fasttree_from_aligned I/O error: %s", exc)
        return "", []
    finally:
        for p in (tmp_fa, tmp_nwk):
            if p and os.path.isfile(p):
                try:
                    os.unlink(p)
                except OSError:
                    pass


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
    """Decompress CDS from cds_seq_map + cds_seq_uniq, trying transcript ID variants."""
    for vid in gene_id_coord_variants(gene_id):
        cur.execute(
            """
            SELECT u.seq_comp
            FROM cds_seq_map m
            JOIN cds_seq_uniq u ON u.uniq_id = m.uniq_id
            WHERE m.gene_id = ? LIMIT 1
            """,
            (vid,),
        )
        r = cur.fetchone()
        if r and r["seq_comp"]:
            try:
                return zlib.decompress(r["seq_comp"]).decode("utf-8")
            except Exception:
                continue
    return None


def gene_cds_map_for_gene_ids(
    conn: sqlite3.Connection, gene_ids: list[str]
) -> dict[str, str]:
    """gene_id -> CDS nucleotide string for genes that have a cds_seq_map row."""
    if not gene_ids or not table_exists(conn, "cds_seq_map"):
        return {}
    cur = conn.cursor()
    out: dict[str, str] = {}
    for gid in gene_ids:
        cds = load_cds_for_gene(cur, gid)
        if cds:
            out[gid] = cds
    return out


def compute_kaks_vs_reference(
    conn,
    pangene_id: str,
    ref_gene_id: str,
    aligned_proteins: dict[str, str],
) -> list[dict[str, Any]]:
    """
    Ka/Ks for each non-reference gene vs ref using CDS mapped onto the given protein alignment.
    aligned_proteins must match the in-app MSA (same gaps as FAMSA output).
    """
    if not table_exists(conn, "cds_seq_map"):
        return []
    ensure_genes_has_pangene(conn)
    cur = conn.cursor()
    cur.execute(f"SELECT gene_id FROM genes WHERE pangene = ?", (pangene_id,))
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
    pangene_id: str,
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
    if not table_exists(conn, "cds_seq_map"):
        return []
    ensure_genes_has_pangene(conn)
    cur = conn.cursor()
    cur.execute(f"SELECT gene_id FROM genes WHERE pangene = ?", (pangene_id,))
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

    Uses batched ``WHERE gene_id IN (...)`` so SQLite can use ``idx_porter6_map_gene_id``.
    (Per-row ``WHERE gene_id = ? COLLATE NOCASE`` cannot use that index and devolves to
    full table scans on large ``porter6_map`` tables.)

    ``porter6_map.gene_id`` may still use ``accession|locus`` in older dumps; rows are
    matched after applying ``canonical_gene_id`` the same way as ``genes.gene_id``.
    """
    out: dict[str, dict[str, str]] = {}
    if not gene_ids:
        return out
    has_table = table_exists(conn, "porter6_map")
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
            expand_seen: set[str] = set()
            expand: list[str] = []
            for g in chunk:
                for x in (g, canonical_gene_id(g)):
                    if x and x not in expand_seen:
                        expand_seen.add(x)
                        expand.append(x)
            ph = ",".join("?" * len(expand))
            cur.execute(
                f"""
                SELECT m.gene_id, u.q3, u.q8
                FROM porter6_map m
                JOIN porter6_uniq u ON u.uniq_id = m.uniq_id
                WHERE m.gene_id IN ({ph})
                """,
                expand,
            )
            for r in cur.fetchall():
                g = (r["gene_id"] or "").strip()
                if not g:
                    continue
                tup = ((r["q3"] or ""), (r["q8"] or ""))
                cg = canonical_gene_id(g)
                qmap[g] = tup
                qmap[cg] = tup
                qmap_lower.setdefault(g.lower(), tup)
                if cg.lower() != g.lower():
                    qmap_lower.setdefault(cg.lower(), tup)

    for gid in gene_ids:
        n = max(1, int(seq_lens.get(gid, 1)))
        if not has_table:
            out[gid] = {"q3": "-" * n, "q8": "-" * n}
            continue
        row = (
            qmap.get(gid)
            or qmap.get(canonical_gene_id(gid))
            or qmap_lower.get(gid.lower())
            or qmap_lower.get(canonical_gene_id(gid).lower())
        )
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


def build_gene_accession_lookup(pangene_id: str, dataset_id: str) -> dict[str, str]:
    """
    Map alignment/FASTA IDs (incl. transcript variants) to genes.accession.
    Wheat IDs are not HORVU-style dotted names; client must not infer accession from split('.')[1].
    """
    m: dict[str, str] = {}
    for g in get_genes_for_pangene(pangene_id, dataset_id):
        acc = g["accession"]
        gid = g["gene_id"]
        for v in gene_id_coord_variants(gid):
            m[v] = acc
    return m


def build_multi_pangene_gene_accession_lookup(
    pangene_ids: list[str], dataset_id: str
) -> dict[str, str]:
    """Union of gene_id → accession for many pan-gene clusters (merged alignment view)."""
    m: dict[str, str] = {}
    conn = get_db(dataset_id)
    try:
        ensure_genes_has_pangene(conn)
        cur = conn.cursor()
        for hid in pangene_ids:
            cur.execute(
                f"SELECT gene_id, accession FROM genes WHERE pangene = ?",
                (hid,),
            )
            for r in cur.fetchall():
                for v in gene_id_coord_variants(r["gene_id"]):
                    m[v] = r["accession"]
    finally:
        conn.close()
    return m


def fetch_gene_row_batch(cur, gene_ids: list[str]) -> dict[str, dict | None]:
    """Gene row per raw ``gene_id``; one SQL query for the whole neighborhood."""
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
    ensure_genes_has_pangene(cur.connection)
    ph = ",".join("?" * len(variants))
    cur.execute(
        f"SELECT gene_id, pangene AS pan, accession FROM genes WHERE gene_id IN ({ph})",
        variants,
    )
    by_norm: dict[str, dict] = {}
    for r in cur.fetchall():
        nk = _norm_gid(r["gene_id"])
        if nk not in by_norm:
            by_norm[nk] = dict(r)
    out: dict[str, dict | None] = {}
    for gid in gene_ids:
        gene_row = None
        for v in gene_id_coord_variants(gid):
            row = by_norm.get(_norm_gid(v))
            if row:
                gene_row = row
                break
        out[gid] = gene_row
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
    chunk_size = 800
    for i in range(0, len(uniq_pairs), chunk_size):
        chunk = uniq_pairs[i : i + chunk_size]
        cur.execute("DROP TABLE IF EXISTS _pv_coord_pairs")
        cur.execute(
            "CREATE TEMP TABLE _pv_coord_pairs "
            "(gene_id TEXT NOT NULL, accession TEXT NOT NULL)"
        )
        cur.executemany(
            "INSERT INTO _pv_coord_pairs (gene_id, accession) VALUES (?, ?)",
            chunk,
        )
        sel_g = ", ".join(f"g.{c.strip()}" for c in sel.split(","))
        cur.execute(
            f"""
            SELECT {sel_g} FROM gene_coords AS g
            INNER JOIN _pv_coord_pairs AS p
              ON g.gene_id = p.gene_id
             AND g.accession = p.accession COLLATE NOCASE
            """
        )
        for r in cur.fetchall():
            rd = dict(r)
            gk = (rd["gene_id"], _norm_acc_key(rd["accession"]))
            fetched[gk] = rd
        cur.execute("DROP TABLE IF EXISTS _pv_coord_pairs")

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


def fetch_gene_row(cur, gene_id: str):
    ensure_genes_has_pangene(cur.connection)
    for vid in gene_id_coord_variants(gene_id):
        cur.execute(
            f"""
            SELECT gene_id, pangene AS pan, accession
            FROM genes WHERE gene_id = ? COLLATE NOCASE
            """,
            (vid,),
        )
        r = cur.fetchone()
        if r:
            return dict(r)
    return None


def _synteny_color_for_key(key: str) -> str:
    """Stable HSL fill for a pan-gene id string.

    Colors do not depend on which other groups appear in the same synteny window,
    so the same pan / type looks identical across stacked rows (e.g. per-gene tracks).
    """
    k = (key or "").strip()
    if not k or k == "—":
        return "#c8c8c8"
    # Discrete palette: more distinct swatches than raw HSL jitter.
    # Chosen to keep decent contrast against the neutral grays used for missing/filtered rows.
    palette = [
        "#1b9e77",
        "#d95f02",
        "#7570b3",
        "#e7298a",
        "#66a61e",
        "#e6ab02",
        "#a6761d",
        "#666666",
        "#1f77b4",
        "#ff7f0e",
        "#2ca02c",
        "#d62728",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#7f7f7f",
        "#bcbd22",
        "#17becf",
        "#393b79",
        "#637939",
        "#843c39",
        "#8c6d31",
        "#ad494a",
        "#f7b6d2",
        "#c7e9c0",
        "#de9f76",
        "#bcbddc",
        "#9ecae1",
        "#fdd0a2",
        "#f768a1",
        "#008080",
        "#6b5b95",
        "#88b04b",
        "#f49ac2",
        "#5b5ea6",
        "#9b2335",
        "#55b4b0",
        "#b565a7",
        "#955251",
        "#009b77",
        "#dd4124",
        "#d65076",
        "#45b8ac",
        "#efc050",
        "#5b4b41",
        "#9b1b30",
        "#009473",
        "#db5640",
        "#743761",
        "#5a7247",
        "#6c4f3d",
        "#587058",
        "#9f9f5c",
        "#b2c248",
        "#407088",
        "#c48d84",
        "#577284",
        "#6b7aa4",
        "#bf9b7b",
        "#4a5366",
        "#7d7f7d",
        "#33658a",
        "#86bbd8",
        "#758ecd",
        "#62466b",
        "#ff6f59",
        "#43aa8b",
        "#f4a259",
        "#bc4b51",
        "#8cb369",
        "#5b9279",
        "#358f80",
        "#e07a5f",
        "#81c14b",
        "#564138",
        "#729ea1",
        "#35524a",
    ]
    digest = hashlib.md5(k.encode("utf-8"), usedforsecurity=False).digest()
    idx = ((digest[0] << 8) | digest[1]) % len(palette)
    return palette[idx]


def build_wheat_cluster_subgenome_table(conn, pangene_id: str, *, limit: int = 600) -> dict:
    """Partition genes in one pan-gene by subgenome (``gene_coords``); may truncate."""
    ensure_genes_has_pangene(conn)
    cur = conn.cursor()
    cur.execute("SELECT COUNT(*) AS c FROM genes WHERE pangene = ?", (pangene_id,))
    total = int(cur.fetchone()["c"])
    cur.execute(
        f"""
        SELECT gene_id, pangene AS pan, accession
        FROM genes WHERE pangene = ?
        ORDER BY accession, gene_id
        LIMIT ?
        """,
        (pangene_id, limit),
    )
    pan_rows = [dict(r) for r in cur.fetchall()]
    by_sg: dict[str, list[dict]] = {"A": [], "B": [], "D": []}
    missing_coords: list[str] = []
    unplaced: list[dict] = []

    for g in pan_rows:
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
        "total_genes": total,
        "shown": len(pan_rows),
        "truncated": total > limit,
    }


def build_wheat_homeologue_table(conn, pangene_id: str) -> dict:
    """A/B/D rows for genes in the same pangene cluster (requires coords for subgenome)."""
    ensure_genes_has_pangene(conn)
    cur = conn.cursor()
    cur.execute(
        f"""
        SELECT gene_id, pangene AS pan, accession
        FROM genes WHERE pangene = ? ORDER BY accession, gene_id
        """,
        (pangene_id,),
    )
    member_genes = [dict(r) for r in cur.fetchall()]
    by_sg: dict[str, list[dict]] = {"A": [], "B": [], "D": []}
    missing_coords: list[str] = []
    unplaced: list[dict] = []

    for g in member_genes:
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


def build_wheat_cluster_pangene_triad_panel(
    conn,
    scope_pangene: str,
    current_pangene_id: str,
    *,
    max_pangenes: int = 40,
) -> dict | None:
    """
    List distinct pan-genes in ``scope_pangene`` (one id → one pan-gene row set), with dominant subgenome.
    When there are exactly 3 pangenes whose dominants are A, B, D, set triad_columns for a 3-column layout.
    """
    if not table_exists(conn, "gene_coords"):
        return None
    ensure_genes_has_pangene(conn)
    cur = conn.cursor()
    cur.execute(
        "SELECT COUNT(DISTINCT pangene) AS n FROM genes WHERE pangene = ?", (scope_pangene,)
    )
    nh = int(cur.fetchone()["n"])
    if nh > max_pangenes:
        return {
            "pangene": scope_pangene,
            "skipped": True,
            "reason": (
                f"This pan-gene scope has {nh} distinct pan-gene rows (showing this panel only when ≤ {max_pangenes})."
            ),
            "pangene_count": nh,
        }

    cur.execute(
        "SELECT DISTINCT pangene AS pan FROM genes WHERE pangene = ? ORDER BY pangene",
        (scope_pangene,),
    )
    pans = [(r["pan"] or "").strip() for r in cur.fetchall()]
    rows: list[dict] = []
    for h in pans:
        cur.execute(
            f"SELECT gene_id, accession FROM genes WHERE pangene = ?",
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
        "pangene": scope_pangene,
        "skipped": False,
        "rows": rows,
        "triad_columns": triad_columns,
    }


def _synteny_strand_is_minus(strand: str | None) -> bool:
    x = (strand or "").strip().lower()
    return x in ("-", "rev", "reverse", "r", "minus")


SYNTENY_NEIGHBOR_WIN_MIN = 1
SYNTENY_FLANK_WIN_MIN = 0
SYNTENY_NEIGHBOR_WIN_MAX = 15


def _clamp_synteny_neighbor_window(n: int) -> int:
    return max(SYNTENY_NEIGHBOR_WIN_MIN, min(int(n), SYNTENY_NEIGHBOR_WIN_MAX))


def _clamp_synteny_flank_window(n: int) -> int:
    """Per-side flank (Collinearity); 0 = focal gene only on that track."""
    return max(SYNTENY_FLANK_WIN_MIN, min(int(n), SYNTENY_NEIGHBOR_WIN_MAX))


def build_wheat_synteny(
    conn,
    focal_gene_id: str,
    focal_accession: str,
    *,
    window: int = 5,
    window_5prime: int | None = None,
    window_3prime: int | None = None,
) -> dict | None:
    """
    Neighboring genes on same chr/accession (genomic order only).
    ``window_5prime`` / ``window_3prime`` are relative to the focal gene strand; if
    both are omitted, symmetric ``window`` is used on each side (2×window + 1 genes).
    Arrow colors use pan-gene id from ``genes.pangene``.
    """
    if window_5prime is None and window_3prime is None:
        w5 = w3 = _clamp_synteny_neighbor_window(window)
    else:
        w5 = _clamp_synteny_flank_window(
            window_5prime if window_5prime is not None else window
        )
        w3 = _clamp_synteny_flank_window(
            window_3prime if window_3prime is not None else window
        )

    cur = conn.cursor()
    fc = fetch_coords_row(cur, focal_gene_id, focal_accession)
    if not fc:
        return None

    chrom = fc["chr"]
    acc = fc["accession"]
    focal_keys = {_norm_gid(x) for x in gene_id_coord_variants(focal_gene_id)}
    minus = _synteny_strand_is_minus(fc.get("strand"))
    w_low = w3 if minus else w5
    w_high = w5 if minus else w3

    slot_count = w_low + w_high + 1
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
        first_ci = max(0, F - w_low)
        last_ci = min(max_ci, F + w_high)
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
            tci = F - w_low + k
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
            j = idx - w_low + k
            slot_js.append(j if 0 <= j < len(chrom_genes) else None)

    if not chrom_genes and not any(j is not None for j in slot_js):
        return None

    present = [chrom_genes[j] for j in slot_js if j is not None]
    if not present:
        return None
    min_s = min(g["start"] for g in present)
    max_e = max(g["end"] for g in present)

    row_batch = fetch_gene_row_batch(
        cur, [g["gene_id"] for g in present]
    )

    missing_rows: list[str] = []
    segments: list[dict] = []
    legend = []
    seen_leg: set[str] = set()

    for j in slot_js:
        if j is None:
            segments.append({"empty": True})
            continue
        g = chrom_genes[j]
        grow = row_batch.get(g["gene_id"])
        if not grow:
            missing_rows.append(g["gene_id"])
            ck = "—"
        else:
            ck = grow["pan"] if grow else "—"
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
                "pangene": grow["pan"] if grow else "",
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
        "color_by": "pan",
        "min_pos": min_s,
        "max_pos": max_e,
        "svg_w": svg_w,
        "svg_h": svg_h,
        "segments": segments,
        "legend": legend,
        "missing_rows": sorted(set(missing_rows)),
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
            f'data-pangene="{escape(str(s.get("pangene") or ""))}">'
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


def enrich_wheat_results_groups(pangene_groups: list[dict], query: str):
    conn = get_db("wheat")
    try:
        if not table_exists(conn, "gene_coords"):
            for g in pangene_groups:
                g["wheat_homeologue"] = None
                g["wheat_cluster_subgenome"] = None
                g["wheat_cluster_pangene_triad"] = None
                g["wheat_synteny"] = None
                g["wheat_synteny_svg"] = ""
            return
        qn = query.strip().lower()

        for group in pangene_groups:
            pan = group["pangene"]
            group["wheat_homeologue"] = build_wheat_homeologue_table(conn, pan)
            group["wheat_cluster_subgenome"] = build_wheat_cluster_subgenome_table(
                conn, pan
            )
            group["wheat_cluster_pangene_triad"] = build_wheat_cluster_pangene_triad_panel(
                conn, pan, pan
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
                    window=5,
                )
            group["wheat_synteny"] = syn
            group["wheat_synteny_svg"] = render_wheat_synteny_svg(syn) if syn else ""
    finally:
        conn.close()


def enrich_wheat_pandagma_results_groups(
    pangene_groups: list[dict],
    query: str,
    *,
    dataset_id: str,
):
    """
    Pandagma mode enrichment:
    - keep only local synteny for the query match window
    - avoid wheat triad panels that assume multi–pan-gene sibling layouts
    """
    conn = get_db(dataset_id)
    try:
        if not table_exists(conn, "gene_coords"):
            for g in pangene_groups:
                g["wheat_homeologue"] = None
                g["wheat_cluster_subgenome"] = None
                g["wheat_cluster_pangene_triad"] = None
                g["wheat_synteny"] = None
                g["wheat_synteny_svg"] = ""
            return

        qn = query.strip().lower()
        for group in pangene_groups:
            group["wheat_homeologue"] = None
            group["wheat_cluster_subgenome"] = None
            group["wheat_cluster_pangene_triad"] = None

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
                    window=5,
                )
            group["wheat_synteny"] = syn
            group["wheat_synteny_svg"] = render_wheat_synteny_svg(syn) if syn else ""
    finally:
        conn.close()


# --- Routes ---

def default_dataset_id() -> str:
    """Prefer wheat when available; otherwise fall back to any installed dataset."""
    d = available_datasets()
    if "wheat" in d:
        return "wheat"
    if d:
        return sorted(d)[0]
    return DEFAULT_DATASET_ID


def _nav_species_ids_ordered() -> list[str]:
    """Wheat, barley, oat (when any variant is installed), then remaining species ids."""
    pref = ["wheat", "barley", "oat"]
    cfg = dataset_configs()
    present = {v["species_id"] for v in cfg.values()}
    out: list[str] = []
    seen: set[str] = set()
    for p in pref:
        if p in present:
            out.append(p)
            seen.add(p)
    for sid in sorted(present):
        if sid not in seen:
            out.append(sid)
    return out


def coerce_dataset_id(raw: str | None = None) -> str:
    """Resolve ``?dataset_id=`` (species tab or variant) to an installed variant id."""
    if raw is None:
        r = (request.args.get("dataset_id") or default_dataset_id()).strip().lower()
    else:
        r = (raw or default_dataset_id()).strip().lower()
    return resolve_dataset_id(r) or default_dataset_id()


def _request_dataset_id() -> str:
    view_args = request.view_args or {}
    if "dataset_id" in view_args:
        path_id = (view_args.get("dataset_id") or "").strip().lower()
        if path_id:
            return coerce_dataset_id(path_id)
    return coerce_dataset_id()


def variant_switch_url(variant_id: str) -> str:
    """Same view with a different index variant (preserves path and query args)."""
    vid = (variant_id or "").strip().lower()
    if not request.endpoint:
        return url_for("dataset", dataset_id=vid)
    view_args = dict(request.view_args or {})
    query = request.args.to_dict(flat=True)
    if "dataset_id" in view_args:
        # ``/dataset/<dataset_id>`` — path param and query must not both pass dataset_id.
        view_args["dataset_id"] = vid
        query.pop("dataset_id", None)
    else:
        query["dataset_id"] = vid
    return url_for(request.endpoint, **view_args, **query)


def _dataset_page_context(
    dataset_id: str, *, show_variant_toggle: bool | None = None
) -> dict[str, Any]:
    cfg = dataset_configs()
    entry = cfg[dataset_id]
    species_id = entry["species_id"]
    variants = variants_for_species(species_id)
    if show_variant_toggle is None:
        show_variant_toggle = request.endpoint in ("index", "search", "dataset")
    return {
        "dataset_id": dataset_id,
        "species_id": species_id,
        "dataset_label": entry["species_label"],
        "variant_label": entry["variant_label"],
        "index_method": entry.get("method") or entry["variant_label"],
        "dataset_variants": variants,
        "about_stats": load_about_stats(dataset_id),
        "show_variant_toggle": show_variant_toggle,
    }


app.jinja_env.globals["variant_switch_url"] = variant_switch_url


@app.context_processor
def inject_nav_defaults():
    dataset_id = _request_dataset_id()
    cfg = dataset_configs()
    species_id = cfg[dataset_id]["species_id"] if dataset_id in cfg else dataset_id
    nav_datasets = [
        {
            "id": sid,
            "label": species_configs()[sid]["label"],
            "default_variant": default_variant_for_species(sid) or sid,
        }
        for sid in _nav_species_ids_ordered()
    ]
    show_toggle = request.endpoint in ("index", "search", "dataset")
    return {
        "nav_default_dataset": default_dataset_id(),
        "nav_datasets": nav_datasets,
        "nav_species_id": species_id,
        "active_dataset_id": dataset_id,
        "dataset_variants": variants_for_species(species_id),
        "variant_label": cfg[dataset_id]["variant_label"] if dataset_id in cfg else "",
        "show_variant_toggle": show_toggle,
        "large_pangene_gene_threshold": LARGE_PANGENE_GENE_THRESHOLD,
    }


@app.route("/dataset/<dataset_id>")
def dataset(dataset_id):
    did = (dataset_id or "").strip().lower()
    resolved = resolve_dataset_id(did)
    if resolved:
        return render_template("index.html", **_dataset_page_context(resolved))
    spec = species_configs().get(did)
    label = spec["label"] if spec else did.capitalize()
    return render_template(
        "dataset_unavailable.html",
        dataset_id=did,
        dataset_label=label,
    )


@app.route("/about")
def about_redirect():
    dataset_id = _request_dataset_id()
    return redirect(url_for("tutorial", dataset_id=dataset_id))


@app.route("/tutorial")
def tutorial():
    dataset_id = _request_dataset_id()
    return render_template("about.html", **_dataset_page_context(dataset_id))


@app.route("/")
def index():
    return redirect(url_for("dataset", dataset_id=default_dataset_id()))


@app.route("/search")
def search():
    dataset_id = coerce_dataset_id()
    species_id = species_id_for_dataset(dataset_id)
    query_raw = request.args.get("q", "").strip()
    query = strip_unsafe_query_chars(query_raw)
    if not query:
        err = (
            "Search query contains unsupported characters."
            if query_raw
            else "Please enter a gene ID."
        )
        return render_template(
            "index.html",
            error=err,
            query=query,
            **_dataset_page_context(dataset_id),
        )

    try:
        results = search_genes(query, dataset_id)
    except (RuntimeError, sqlite3.OperationalError) as e:
        logging.getLogger(__name__).exception(
            "gene search failed q=%r dataset_id=%s", query, dataset_id
        )
        page_ctx = _dataset_page_context(dataset_id, show_variant_toggle=True)
        return render_template(
            "results.html",
            query=query,
            error=f'Search failed for "{query}".',
            dataset_id=dataset_id,
            species_id=species_id,
            index_method=page_ctx.get("index_method"),
            pangene_groups=[],
            show_wheat_extra=False,
            wheat_coords_available=False,
            coords_available=False,
            keyword_mode=False,
            matched_refs=[],
            keyword_selected_refs=[],
        )
    # Wrong dataset is the usual cause: default is wheat when both DBs exist, but ``HORVU.`` hits barley.
    if not results and query:
        barley_v = default_variant_for_species("barley")
        wheat_v = default_variant_for_species("wheat")
        oat_v = default_variant_for_species("oat")
        if (
            species_id == "wheat"
            and _looks_like_barley_gene_query(query)
            and barley_v
        ):
            if search_genes(query, barley_v):
                return redirect(
                    url_for("search", q=query, dataset_id=barley_v)
                )
        if (
            species_id == "barley"
            and _looks_like_wheat_gene_query(query)
            and wheat_v
        ):
            if search_genes(query, wheat_v):
                return redirect(url_for("search", q=query, dataset_id=wheat_v))
        if (
            species_id in ("wheat", "barley")
            and _looks_like_oat_gene_query(query)
            and oat_v
        ):
            if search_genes(query, oat_v):
                return redirect(url_for("search", q=query, dataset_id=oat_v))
        if (
            species_id == "oat"
            and _looks_like_wheat_gene_query(query)
            and wheat_v
        ):
            if search_genes(query, wheat_v):
                return redirect(url_for("search", q=query, dataset_id=wheat_v))
        if species_id == "oat" and _looks_like_barley_gene_query(query) and barley_v:
            if search_genes(query, barley_v):
                return redirect(url_for("search", q=query, dataset_id=barley_v))

    # Keyword fallback: when the query does not resolve to any pan-genome gene
    # ID (after cross-species redirects), try the shared search index built
    # from Arabidopsis/rice annotations joined via mmseqs best hits.
    keyword_refs: list[dict] = []
    keyword_pangene_groups: list[dict] = []
    keyword_selected_refs: list[str] = []
    if not results and query and keyword_index_available():
        _cfg_kw = get_dataset_config(dataset_id)
        # Accept repeated ``?refs=A&refs=B`` (checkbox form) or legacy
        # comma-separated ``?refs=A,B`` (bookmarkable URL).
        refs_raw = request.args.getlist("refs")
        if len(refs_raw) == 1 and "," in refs_raw[0]:
            refs_raw = refs_raw[0].split(",")
        keyword_selected_refs = [x.strip() for x in refs_raw if x and x.strip()]
        restrict_refs: set[str] | None = (
            set(keyword_selected_refs) if keyword_selected_refs else None
        )
        try:
            kw = search_keywords(
                query,
                keyword_dataset_id_for_dataset(dataset_id),
                db_path=_cfg_kw.get("db_path"),
                restrict_refs=restrict_refs,
            )
        except RuntimeError as e:
            logging.getLogger(__name__).exception(
                "keyword search failed q=%r dataset_id=%s", query, dataset_id
            )
            page_ctx = _dataset_page_context(dataset_id, show_variant_toggle=True)
            return render_template(
                "results.html",
                query=query,
                error=f'Search failed for "{query}".',
                dataset_id=dataset_id,
                species_id=species_id,
                index_method=page_ctx.get("index_method"),
                pangene_groups=[],
                show_wheat_extra=False,
                wheat_coords_available=False,
                coords_available=False,
                keyword_mode=False,
                matched_refs=[],
                keyword_selected_refs=[],
            )
        except (sqlite3.Error, OSError, ValueError, TypeError):
            logging.getLogger(__name__).exception(
                "keyword search failed q=%r dataset_id=%s refs=%s",
                query,
                dataset_id,
                keyword_selected_refs,
            )
            page_ctx = _dataset_page_context(dataset_id, show_variant_toggle=True)
            return render_template(
                "results.html",
                query=query,
                error=f'Search failed for "{query}".',
                dataset_id=dataset_id,
                species_id=species_id,
                index_method=page_ctx.get("index_method"),
                pangene_groups=[],
                show_wheat_extra=False,
                wheat_coords_available=False,
                coords_available=False,
                keyword_mode=False,
                matched_refs=[],
                keyword_selected_refs=keyword_selected_refs,
            )
        keyword_refs = kw.get("refs", [])
        keyword_pangene_groups = kw.get("pangene_groups", [])

    if not results and not keyword_refs and not keyword_pangene_groups:
        page_ctx = _dataset_page_context(dataset_id, show_variant_toggle=True)
        return render_template(
            "results.html",
            query=query,
            error=f'No results found for "{query}".',
            dataset_id=dataset_id,
            species_id=species_id,
            index_method=page_ctx.get("index_method"),
            pangene_groups=[],
            show_wheat_extra=False,
            wheat_coords_available=False,
            coords_available=False,
            keyword_mode=False,
            matched_refs=[],
            keyword_selected_refs=[],
        )

    # Keyword mode is on whenever the keyword index returned anything for the
    # current dataset, even before the user has selected reference genes.
    keyword_mode = (not results) and (
        bool(keyword_refs) or bool(keyword_pangene_groups)
    )

    if results and len(results) == 1:
        gene = results[0]
        n_genes = count_pangene_genes(gene["pan"], dataset_id)
        if n_genes < LARGE_PANGENE_GENE_THRESHOLD or request.args.get("ack_large"):
            detail_qs: dict[str, str] = {
                "highlight": query,
                "dataset_id": dataset_id,
            }
            if request.args.get("ack_large"):
                detail_qs["ack_large"] = "1"
            return redirect(
                url_for(
                    "pangene_detail",
                    pangene_id=gene["pan"],
                    **detail_qs,
                )
            )

    if keyword_mode:
        pangene_list = keyword_pangene_groups
    else:
        pans_seen: dict[str, dict] = {}
        for r in results:
            pan = r["pan"]
            if pan not in pans_seen:
                pans_seen[pan] = {
                    "pangene": pan,
                    "genes": [],
                    "gene_count": count_pangene_genes(pan, dataset_id),
                }
            pans_seen[pan]["genes"].append(
                {"gene_id": r["gene_id"], "accession": r["accession"]}
            )
        pangene_list = list(pans_seen.values())

    coords_available = False
    try:
        _conn = get_db(dataset_id)
    except (RuntimeError, sqlite3.OperationalError) as e:
        logging.getLogger(__name__).warning(
            "search coords skipped (db open failed) dataset_id=%s: %s",
            dataset_id,
            e,
        )
        _conn = None
    if _conn is not None:
        try:
            coords_available = table_exists(_conn, "gene_coords")
            if keyword_mode and coords_available and pangene_list:
                cur = _conn.cursor()
                ga_pairs: list[tuple[str, str]] = []
                for grp in pangene_list:
                    for ge in grp.get("genes", []):
                        ga_pairs.append((ge["gene_id"], ge["accession"]))
                if ga_pairs:
                    row_map = fetch_coords_rows_batch(cur, ga_pairs)
                    for grp in pangene_list:
                        for ge in grp.get("genes", []):
                            cr = row_map.get(ge["gene_id"])
                            if cr:
                                ge["chr"] = cr.get("chr")
                                ge["start"] = cr.get("start")
                                ge["end"] = cr.get("end")
                                ge["strand"] = cr.get("strand")
        except (RuntimeError, sqlite3.Error, OSError) as e:
            logging.getLogger(__name__).warning(
                "search coords lookup skipped q=%r dataset_id=%s: %s",
                query,
                dataset_id,
                e,
            )
            coords_available = False
        finally:
            _conn.close()

    if coords_available and not keyword_mode:
        cfg = get_dataset_config(dataset_id)
        if cfg.get("mode") == "pandagma" and not cfg.get("enable_homeologue_panels"):
            enrich_wheat_pandagma_results_groups(
                pangene_list, query, dataset_id=dataset_id
            )
        elif species_id == "wheat":
            enrich_wheat_results_groups(pangene_list, query)

    page_ctx = _dataset_page_context(dataset_id, show_variant_toggle=True)
    return render_template(
        "results.html",
        query=query,
        dataset_id=dataset_id,
        species_id=species_id,
        index_method=page_ctx.get("index_method"),
        pangene_groups=pangene_list,
        show_wheat_extra=(species_id == "wheat"),
        wheat_coords_available=coords_available,
        coords_available=coords_available,
        keyword_mode=keyword_mode,
        matched_refs=keyword_refs if keyword_mode else [],
        keyword_selected_refs=keyword_selected_refs if keyword_mode else [],
    )


@app.route("/pangene/<pangene_id>")
def pangene_detail(pangene_id):
    dataset_id = coerce_dataset_id()
    species_id = species_id_for_dataset(dataset_id)
    highlight = request.args.get("highlight", "")
    pan_internal = internal_pangene_id(pangene_id, dataset_id)
    info = get_pangene_info(pan_internal, dataset_id)
    if not info:
        return render_template(
            "index.html",
            error=f'Pangene "{pangene_id}" not found.',
            **_dataset_page_context(dataset_id),
        )

    pangene_id = pan_internal
    gene_count = count_pangene_genes(pangene_id, dataset_id)
    if gene_count >= LARGE_PANGENE_GENE_THRESHOLD and not request.args.get("ack_large"):
        proceed_qs = {"dataset_id": dataset_id, "ack_large": "1"}
        if highlight:
            proceed_qs["highlight"] = highlight
        return render_template(
            "large_pangene_confirm.html",
            pangene_id=pangene_id,
            pangene_label=pangene_display_label(pangene_id, dataset_id),
            gene_count=gene_count,
            proceed_url=url_for("pangene_detail", pangene_id=pangene_id, **proceed_qs),
            **_dataset_page_context(dataset_id),
        )

    genes = get_genes_for_pangene(pangene_id, dataset_id)
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

    wheat_cluster_picker = None
    wheat_synteny_available = False
    gene_row_meta: dict[str, dict] = {}
    pandagma_pan_view_available = False
    cds_available = False
    porter6_by_gene: dict[str, dict[str, str]] = {}
    _conn = get_db(dataset_id)
    try:
        cds_available = table_exists(_conn, "cds_seq_map")
        wheat_synteny_available = table_exists(_conn, "gene_coords")
        cfg = get_dataset_config(dataset_id)
        pandagma_pan_view_available = (
            cfg.get("mode") == "pandagma" and wheat_synteny_available
        )
        if dataset_id == "wheat" and get_dataset_config(dataset_id).get(
            "enable_cluster_picker"
        ):
            wheat_cluster_picker = build_wheat_cluster_picker_payload(
                _conn, pangene_id, pangene_id
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
                    m["subgenome"] = subgenome_for_ui_from_coords_row(cr, dataset_id)
            gene_row_meta[gid] = m
        porter6_by_gene = (
            porter6_raw_for_aligned_sequences(_conn, aligned) if aligned else {}
        )
    finally:
        _conn.close()

    gene_seq_lookup: dict[str, dict[str, str]] = {"cds": {}, "protein": {}}
    if aligned:
        for gid in aligned:
            raw_p = sequences.get(gid)
            if raw_p:
                gene_seq_lookup["protein"][gid] = raw_p
        if cds_available:
            _conn2 = get_db(dataset_id)
            try:
                gene_seq_lookup["cds"] = gene_cds_map_for_gene_ids(
                    _conn2, list(aligned.keys())
                )
            finally:
                _conn2.close()

    chinese_spring_gene_ids = chinese_spring_gene_ids_from_rows(genes)
    morex_v3_gene_ids = morex_v3_reference_gene_ids_from_rows(genes)
    sang_v11_gene_ids = sang_v11_reference_gene_ids_from_rows(genes)
    if species_id == "wheat":
        expression_panel_enabled = True
        expression_reference_gene_ids = chinese_spring_gene_ids
        expression_ref_accession = "chinesespring"
        expression_plantapp_genome = ""
        expression_ref_label = "Chinese Spring"
    elif species_id == "barley":
        expression_panel_enabled = True
        expression_reference_gene_ids = morex_v3_gene_ids
        expression_ref_accession = "morexv3"
        expression_plantapp_genome = "HvMorex"
        expression_ref_label = "MorexV3 (PlantApp: HvMorex)"
    elif species_id == "oat":
        expression_panel_enabled = True
        expression_reference_gene_ids = sang_v11_gene_ids
        expression_ref_accession = "sangv11"
        expression_plantapp_genome = "AsSang"
        expression_ref_label = "sangV11 (PlantApp: AsSang)"
    else:
        expression_panel_enabled = False
        expression_reference_gene_ids = []
        expression_ref_accession = ""
        expression_plantapp_genome = ""
        expression_ref_label = ""

    register_latest_famsa_alignment(aligned)

    gene_ids_for_cross: list[str] = list(
        dict.fromkeys(
            [
                str(row["gene_id"])
                for row in genes
                if row["gene_id"] is not None and str(row["gene_id"]).strip()
            ]
        )
    )
    if aligned:
        for gid in aligned:
            if gid not in gene_ids_for_cross:
                gene_ids_for_cross.append(gid)
    cross_species_by_gene: dict[str, Any] = {}
    if keyword_index_available():
        try:
            cross_species_by_gene = cross_species_best_annotations_for_genes(
                keyword_dataset_id_for_dataset(dataset_id), gene_ids_for_cross
            )
        except (RuntimeError, OSError, sqlite3.Error, TypeError, ValueError):
            cross_species_by_gene = {}

    return render_template(
        "pangene_detail.html",
        dataset_id=dataset_id,
        pangene_id=pangene_id,
        info=info,
        genes=genes,
        gene_row_meta=gene_row_meta,
        sequences=sequences,
        gene_seq_lookup=gene_seq_lookup,
        aligned_json=json.dumps(aligned),
        consensus_str=consensus_str,
        aln_len=aln_len,
        highlight=highlight,
        num_seqs=len(aligned),
        newick=newick,
        leaf_order_json=json.dumps(leaf_order),
        conservation_json=json.dumps(conservation_scores),
        wheat_cluster_picker=wheat_cluster_picker,
        wheat_synteny_available=wheat_synteny_available,
        gene_accession_lookup=gene_accession_lookup,
        pandagma_pan_view_available=pandagma_pan_view_available,
        cds_available=cds_available,
        porter6_json=porter6_by_gene,
        external_db_links=external_db_links_for_dataset(dataset_id),
        genome_metadata=genome_metadata_for_dataset(dataset_id),
        grain_genes_blast=grain_genes_blast_config_for_dataset(dataset_id),
        chinese_spring_gene_ids=chinese_spring_gene_ids,
        expression_panel_enabled=expression_panel_enabled,
        expression_reference_gene_ids=expression_reference_gene_ids,
        expression_ref_accession=expression_ref_accession,
        expression_plantapp_genome=expression_plantapp_genome,
        expression_ref_label=expression_ref_label,
        genome_subgenome_layout=genome_subgenome_layout(dataset_id),
        gene_subgenome_toggle_label=(
            "A / B / D"
            if genome_subgenome_layout(dataset_id) == "abd"
            else (
                "A / C / D"
                if genome_subgenome_layout(dataset_id) == "acd"
                else "Chr"
            )
        ),
        cross_species_by_gene=cross_species_by_gene,
    )


@app.route("/api/plantapp_omics")
def api_plantapp_omics():
    """
    Proxy PlantApp tissue + DEG JSON for one gene_id (see PlantApp pages/api.py).
    Optional ``genome`` (e.g. ``HvMorex`` for barley) is forwarded when provided.
    """
    gene_id = (request.args.get("gene_id") or "").strip()
    genome = (request.args.get("genome") or "").strip() or None
    if not gene_id:
        return jsonify({"error": "gene_id is required"}), 400
    try:
        return jsonify(fetch_plantapp_omics(gene_id, genome=genome))
    except Exception:
        logging.getLogger(__name__).exception(
            "plantapp_omics failed gene_id=%r genome=%r", gene_id, genome
        )
        return (
            jsonify(
                {
                    "query_gene_id": gene_id,
                    "ok": False,
                    "error": "PlantApp expression request failed on the server.",
                }
            ),
            500,
        )


@app.route("/pangene/<pangene_id>/tree")
def pangene_tree(pangene_id):
    dataset_id = coerce_dataset_id()
    pan_internal = internal_pangene_id(pangene_id, dataset_id)
    info = get_pangene_info(pan_internal, dataset_id)
    if not info:
        return render_template(
            "index.html",
            error=f'Pangene "{pangene_id}" not found.',
            **_dataset_page_context(dataset_id),
        )

    pangene_id = pan_internal
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

    register_latest_famsa_alignment(aligned)

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
    dataset_id = coerce_dataset_id()
    export_stamp = date.today().strftime("%Y%m%d")

    def export_filename(suffix: str) -> str:
        safe_id = re.sub(r"[^\w\-.]+", "_", str(pangene_id or "pangene"))
        return f"{safe_id}.panviewer.{export_stamp}.{suffix}"

    info = get_pangene_info(pangene_id, dataset_id)
    if not info:
        return "Pangene not found", 404

    genes_all = get_genes_for_pangene(pangene_id, dataset_id)
    cluster_ids = {str(g["gene_id"]) for g in genes_all}
    requested = [x.strip() for x in request.args.getlist("gene_id") if x.strip()]
    allowed_order: list[str] | None = None
    if requested:
        allowed_order = [g for g in requested if g in cluster_ids]
        if not allowed_order:
            return (
                "None of the requested gene IDs belong to this pan-gene cluster",
                400,
            )

    if dtype == "gene_list":
        genes = genes_all
        if allowed_order is not None:
            by_gid = {str(g["gene_id"]): g for g in genes_all}
            genes = [by_gid[g] for g in allowed_order if g in by_gid]
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
                "Content-Disposition": f"attachment; filename={export_filename('genes.tsv')}"
            },
        )

    elif dtype in ("proteins", "sequences"):
        sequences = read_fasta(pangene_id, dataset_id)
        if allowed_order is not None:
            sequences = {
                g: sequences[g] for g in allowed_order if g in sequences
            }
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
                "Content-Disposition": f"attachment; filename={export_filename('proteins.fasta')}"
            },
        )

    elif dtype == "cds":
        conn = get_db(dataset_id)
        try:
            if not table_exists(conn, "cds_seq_map"):
                return "CDS sequences are not available for this dataset", 503
            gene_ids = (
                allowed_order
                if allowed_order is not None
                else [str(g["gene_id"]) for g in genes_all]
            )
            cds_map = gene_cds_map_for_gene_ids(conn, gene_ids)
        finally:
            conn.close()
        lines = []
        for gid in gene_ids:
            seq = cds_map.get(gid)
            if not seq:
                continue
            lines.append(f">{gid}")
            for i in range(0, len(seq), 80):
                lines.append(seq[i : i + 80])
        content = "\n".join(lines)
        return Response(
            content,
            mimetype="text/plain",
            headers={
                "Content-Disposition": f"attachment; filename={export_filename('cds.fasta')}"
            },
        )

    elif dtype == "alignment":
        sequences = read_fasta(pangene_id, dataset_id)
        if allowed_order is not None:
            sequences = {
                g: sequences[g] for g in allowed_order if g in sequences
            }
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
                "Content-Disposition": f"attachment; filename={export_filename('alignment.fasta')}"
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
                "Content-Disposition": f"attachment; filename={export_filename('variants.tsv')}"
            },
        )

    return "Unknown download type", 400


@app.route("/api/search_suggestions")
def search_suggestions():
    dataset_id = coerce_dataset_id()
    query = strip_unsafe_query_chars(request.args.get("q", "").strip())
    if len(query) < 3:
        return jsonify([])
    like_pat = f"%{escape_sql_like(canonical_gene_id(query))}%"
    conn = get_db(dataset_id)
    ensure_genes_has_pangene(conn)
    cur = conn.cursor()
    meta = pangene_info_meta(conn)
    if meta:
        tbl, col = meta
        cur.execute(
            f"SELECT g.gene_id, g.pangene AS pan, g.accession, "
            f"COALESCE(pi.gene_count, (SELECT COUNT(*) FROM genes g2 WHERE g2.pangene = g.pangene)) "
            f"AS gene_count FROM genes g "
            f"LEFT JOIN {tbl} pi ON pi.{col} = g.pangene "
            "WHERE g.gene_id LIKE ? ESCAPE '\\' COLLATE NOCASE LIMIT 15",
            (like_pat,),
        )
    else:
        cur.execute(
            f"SELECT gene_id, pangene AS pan, accession, "
            "(SELECT COUNT(*) FROM genes g2 WHERE g2.pangene = genes.pangene) AS gene_count "
            "FROM genes WHERE gene_id LIKE ? ESCAPE '\\' COLLATE NOCASE LIMIT 15",
            (like_pat,),
        )
    results = [
        {
            "gene_id": r["gene_id"],
            "display_gene_id": gene_id_display_label(r["gene_id"]),
            "pangene": r["pan"],
            "accession": r["accession"],
            "gene_count": int(r["gene_count"] or 0),
        }
        for r in cur.fetchall()
    ]
    conn.close()
    return jsonify(results)


@app.route("/api/wheat/synteny_tracks", methods=["POST"])
def api_wheat_synteny_tracks():
    """One ±window synteny track per pangene (focal = first gene with coords in that cluster)."""
    data = request.get_json(silent=True) or {}
    dataset_id = coerce_dataset_id(data.get("dataset_id"))
    if dataset_id not in available_datasets():
        return jsonify(error="dataset not available"), 404
    conn = get_db(dataset_id)
    try:
        if not table_exists(conn, "gene_coords"):
            return jsonify(error="gene_coords not loaded"), 503
        raw_ids = data.get("pangene_ids") or []
        if not isinstance(raw_ids, list):
            return jsonify(error="pangene_ids must be a list"), 400
        pangene_ids = [str(h).strip() for h in raw_ids if str(h).strip()]
        cur = conn.cursor()
        ensure_genes_has_pangene(conn)
        tracks: list[dict] = []
        for hid in pangene_ids:
            cur.execute(
                f"SELECT gene_id, accession FROM genes WHERE pangene = ? ORDER BY gene_id",
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
                    "missing_rows": syn.get("missing_rows", []),
                }
            )
        return jsonify(tracks=tracks)
    finally:
        conn.close()


@app.route("/api/wheat/merged_alignment", methods=["POST"])
def api_wheat_merged_alignment():
    """FAMSA (NJ guide tree) across proteins for the listed pan-gene clusters (wheat)."""
    if "wheat" not in available_datasets():
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
    acc_map = build_multi_pangene_gene_accession_lookup(pangene_ids, "wheat")

    conn = get_db("wheat")
    try:
        porter6 = porter6_raw_for_aligned_sequences(conn, aligned)
        gene_cds = gene_cds_map_for_gene_ids(conn, list(aligned.keys()))
    finally:
        conn.close()

    register_latest_famsa_alignment(aligned)

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
        gene_cds=gene_cds,
    )


@app.route("/api/grain_genes_blast_transfer", methods=["POST"])
def api_grain_genes_blast_transfer():
    """Proxy GrainGenes BLAST ``query_transfer`` (avoids browser CORS on local dev)."""
    payload = request.get_json(silent=True) or {}
    database = (payload.get("database") or "").strip()
    query = payload.get("query")
    if not database or query is None or not str(query).strip():
        return jsonify({"error": "database and query are required"}), 400
    body = json.dumps({"database": database, "query": str(query)}).encode("utf-8")
    req = urllib.request.Request(
        GRAINGENES_BLAST_TRANSFER_URL,
        data=body,
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            return jsonify(json.loads(resp.read().decode("utf-8")))
    except urllib.error.HTTPError as e:
        detail = e.read().decode("utf-8", errors="replace")[:500]
        return jsonify({"error": f"GrainGenes BLAST transfer failed: HTTP {e.code}", "detail": detail}), 502
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as e:
        return jsonify({"error": f"GrainGenes BLAST transfer failed: {e}"}), 502


@app.route("/api/align_genes", methods=["POST"])
def api_align_genes():
    """Re-align only the requested gene IDs (active subset) with FAMSA."""
    data = request.get_json(silent=True) or {}
    dataset_id = coerce_dataset_id(data.get("dataset_id"))
    if dataset_id not in available_datasets():
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
        gene_cds = gene_cds_map_for_gene_ids(conn2, list(aligned.keys()))
    finally:
        conn2.close()

    register_latest_famsa_alignment(aligned)

    return jsonify(
        aligned=aligned,
        aln_len=aln_len,
        newick=newick or "",
        leaf_order=lo,
        conservation=cons,
        gene_accession=out_acc,
        porter6=porter6,
        num_seqs=len(aligned),
        gene_cds=gene_cds,
    )


@app.route("/api/newick_fasttree", methods=["POST"])
def api_newick_fasttree():
    """FastTree + midpoint on session-cached MSA, or rebuild FAMSA if the cache missed."""
    data = request.get_json(silent=True) or {}
    use_sess = bool(data.get("use_session_msa"))
    aligned: dict[str, str] = {}
    if use_sess:
        dataset_id = coerce_dataset_id(data.get("dataset_id"))
        if dataset_id not in available_datasets():
            return jsonify(error="dataset_not_available"), 404
        raw_ids = data.get("gene_ids")
        if not isinstance(raw_ids, list):
            return jsonify(error="gene_ids_required"), 400
        seen: set[str] = set()
        gene_ids_clean: list[str] = []
        for g in raw_ids:
            s = str(g).strip()
            if s and s not in seen:
                seen.add(s)
                gene_ids_clean.append(s)
        if len(gene_ids_clean) < 2:
            return jsonify(error="need_two_gene_ids"), 400

        hit = get_latest_famsa_alignment_from_session()
        if hit and len(hit) >= 2:
            for g in gene_ids_clean:
                if g in hit:
                    aligned[g] = hit[g]

        if len(aligned) < 2:
            full = famsa_rebuild_and_cache_for_genes(dataset_id, gene_ids_clean)
            if not full:
                return jsonify(error="realign_failed"), 502
            aligned = {g: full[g] for g in gene_ids_clean if g in full}
            if len(aligned) < 2:
                return jsonify(error="realign_missing_sequences"), 400
    else:
        aligned = _sanitize_aligned_for_fasttree_api(data.get("aligned"))
        if len(aligned) < 2:
            return jsonify(error="need_at_least_two_aligned_sequences"), 400
    if len(aligned) > 600:
        return jsonify(error="too_many_sequences_max_600"), 400
    newick, lo = newick_fasttree_from_aligned(aligned)
    if not newick:
        return jsonify(error="fasttree_failed", detail="FastTree missing or could not read the MSA"), 502
    return jsonify(ok=True, newick=newick, leaf_order=lo)


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


@app.route("/api/cross_species_annotations", methods=["POST"])
def api_cross_species_annotations():
    """Best Arabidopsis/rice keyword-index hit per gene_id (for synteny tooltips)."""
    data = request.get_json(silent=True) or {}
    dataset_id = coerce_dataset_id(data.get("dataset_id"))
    if dataset_id not in available_datasets():
        return jsonify(error="dataset not available"), 404
    gene_ids_in = data.get("gene_ids") or []
    if not isinstance(gene_ids_in, list):
        return jsonify(error="gene_ids must be a list"), 400
    gene_ids = [str(g).strip() for g in gene_ids_in if str(g).strip()][:2000]
    if not gene_ids:
        return jsonify({})
    if not keyword_index_available():
        return jsonify({})
    try:
        ann = cross_species_best_annotations_for_genes(
            keyword_dataset_id_for_dataset(dataset_id), gene_ids
        )
        return jsonify(ann)
    except (RuntimeError, OSError, sqlite3.Error, TypeError, ValueError) as e:
        return jsonify(error=str(e)), 500


@app.route("/api/wheat/synteny_per_gene", methods=["GET", "POST", "OPTIONS"])
def api_wheat_synteny_per_gene():
    """One compact equal-width synteny strip per (gene_id, accession) focal."""
    if request.method == "OPTIONS":
        r = Response(status=204)
        r.headers["Allow"] = "GET, POST, OPTIONS"
        return r
    if request.method == "GET":
        return jsonify(
            message="This endpoint requires POST with JSON: genes, dataset_id.",
        )
    data = request.get_json(silent=True) or {}
    dataset_id = coerce_dataset_id(data.get("dataset_id"))
    if dataset_id not in available_datasets():
        return jsonify(error="dataset not available"), 404
    conn = get_db(dataset_id)
    try:
        if not table_exists(conn, "gene_coords"):
            return jsonify(error="gene_coords not loaded"), 503
        genes_in = data.get("genes") or []
        if not isinstance(genes_in, list):
            return jsonify(error="genes must be a list"), 400
        if len(genes_in) > 150:
            return jsonify(error="too_many_genes_max_150"), 400

        try:
            window = int(data.get("window", 5))
        except (TypeError, ValueError):
            window = 5
        window = _clamp_synteny_neighbor_window(window)

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
            w5_item = item.get("window_5prime")
            w3_item = item.get("window_3prime")
            if w5_item is not None or w3_item is not None:
                try:
                    w5v = _clamp_synteny_flank_window(
                        w5_item if w5_item is not None else window
                    )
                except (TypeError, ValueError):
                    w5v = window
                try:
                    w3v = _clamp_synteny_flank_window(
                        w3_item if w3_item is not None else window
                    )
                except (TypeError, ValueError):
                    w3v = window
                syn = build_wheat_synteny(
                    conn, gid, acc, window=window, window_5prime=w5v, window_3prime=w3v
                )
            else:
                syn = build_wheat_synteny(conn, gid, acc, window=window)
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
                    "missing_rows": syn.get("missing_rows", []),
                }
            )

        return Response(json.dumps(tracks), mimetype="application/json")
    finally:
        conn.close()


@app.route("/api/micro_coll_track_order", methods=["POST", "OPTIONS"])
def api_micro_coll_track_order():
    """Order Collinearity rows by protein k-mer similarity (UPGMA on 1 − Jaccard)."""
    if request.method == "OPTIONS":
        r = Response(status=204)
        r.headers["Allow"] = "POST, OPTIONS"
        return r
    data = request.get_json(silent=True) or {}
    genes_in = data.get("genes") or []
    seq_in = data.get("sequences") or {}
    if not isinstance(genes_in, list) or not isinstance(seq_in, dict):
        return jsonify(error="genes_sequences_required"), 400
    if len(genes_in) < 2 or len(genes_in) > 150:
        return jsonify(error="gene_count_must_be_2_to_150"), 400
    gene_ids: list[str] = []
    for g in genes_in:
        if g is None:
            continue
        s = str(g).strip()
        if not s:
            return jsonify(error="empty_gene_id"), 400
        gene_ids.append(s)
    if len(gene_ids) != len(genes_in):
        return jsonify(error="invalid_gene_entries"), 400

    sequences: dict[str, str] = {}
    total_chars = 0
    for gid in gene_ids:
        raw = seq_in.get(gid)
        if not isinstance(raw, str):
            return jsonify(error="sequence_must_be_string", gene_id=gid), 400
        seq = "".join(c for c in raw.strip().upper() if c.isalpha() or c == "-")
        sequences[gid] = seq
        total_chars += len(seq)
    if gene_ids:
        first_len = len(sequences[gene_ids[0]])
        if first_len > 15_000 or total_chars > 800_000:
            return jsonify(error="sequences_too_large"), 400

    order = protein_track_order_permutation(gene_ids, sequences)
    return jsonify(order=order, n=len(order))


@app.route("/api/pairwise_kaks", methods=["POST"])
def api_pairwise_kaks():
    """
    Pairwise Ka/Ks vs a reference from CDS mapped onto the client-provided protein alignment.
    """
    data = request.get_json(silent=True) or {}
    dataset_id = coerce_dataset_id(data.get("dataset_id"))
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
        if not table_exists(conn, "cds_seq_map"):
            return jsonify(error="cds_not_loaded"), 503
        cur = conn.cursor()
        ensure_genes_has_pangene(conn)
        cur.execute(f"SELECT gene_id FROM genes WHERE pangene = ?", (pangene_id,))
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
    if not dataset_configs():
        print(
            "No SQLite files under database/. Run `python build_index.py` "
            "(or `python build_index.py --force`) after placing species inputs under input/."
        )
        exit(1)
    _dd = default_dataset_id()
    default_db = dataset_configs()[_dd]["db_path"]
    if not os.path.exists(default_db):
        print("Default dataset database not found. Run `python build_index.py --force`.")
        exit(1)
    app.run(debug=True, port=5050)
