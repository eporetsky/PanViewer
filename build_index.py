"""
Build per-species SQLite indexes from Pandagma pan-gene membership + BED coordinates,
deduplicated embedded protein/CDS FASTA (proteins: all ``*`` removed), per-cluster metadata
+ gene rows in SQLite, and deduplicated porter6 secondary structure predictions.
"""
import csv
import functools
import gzip
import io
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import zlib

from Bio import SeqIO

from dataset_stats import merge_write_combined_stats_tsv, stats_from_cursor
from gene_id_normalize import canonical_gene_id, fasta_header_token_candidates


def canonical_pan_id_for_db(raw_pan_id: str, species_id: str) -> str:
    """
    Map Pandagma pan column → stored ``genes.pangene`` / ``pangene_info.pangene`` id.

    Wheat: ``pan00001`` → ``Traes_pan00001``. Barley: ``pan00001`` → ``HORVU_pan00001``.
    Oat: ``pan_00001`` → ``Avena_pan00001`` (same ``pan_`` slice rule as wheat UI).
    Other species: unchanged. Skips re-prefixing if the id already carries the species prefix.
    """
    pan = (raw_pan_id or "").strip()
    if not pan:
        return pan
    sid = (species_id or "").strip().lower()
    if sid == "wheat":
        if pan.lower().startswith("traes_"):
            return pan
        return f"Traes_{pan}"
    if sid == "barley":
        if pan.upper().startswith("HORVU_"):
            return pan
        return f"HORVU_{pan}"
    if sid == "oat":
        pl = pan.lower()
        if pl.startswith("avena_"):
            return pan
        if len(pan) >= 4 and pl.startswith("pan_"):
            return "Avena_pan" + pan[4:]
        return f"Avena_{pan}"
    return pan


def subgenome_from_chromosome(chrom: str) -> str | None:
    """Infer subgenome letter from chromosome names (e.g. wheat 1A/7B; oat hexaploid 1A/1C/4D)."""
    if not chrom:
        return None
    s = str(chrom).strip().upper()
    # BED sometimes uses "chr1A" / "chr4D" naming; coords TSV often uses "1A" / "4D".
    if s.startswith("CHR"):
        s = s[3:]
    # Supported patterns:
    # - "1A" / "4D" / "1C"  -> letter at the end
    # - "A1" / "D7" / "C3"  -> letter at the start (e.g., chrA1)
    m = re.fullmatch(r"\d+([ABCD])", s)
    if m:
        return m.group(1)
    m = re.fullmatch(r"([ABCD])\d+", s)
    if m:
        return m.group(1)
    return None


def resolve_fasta_record_to_db_gene_id(
    accession: str, record, valid_gene_ids: set[str]
) -> str | None:
    """Map a FASTA record header to ``genes.gene_id`` (canonical ids only in ``genes``)."""
    _ = accession
    for tok in fasta_header_token_candidates(record):
        cand = canonical_gene_id(tok)
        if cand and cand in valid_gene_ids:
            return cand
    return None


# Longest first so ``.cds.fa.gz`` beats ``.fa.gz``.
_FASTA_FILENAME_SUFFIXES = (
    ".cds.fa.gz",
    ".pep.fa.gz",
    ".faa.gz",
    ".fasta.gz",
    ".fna.gz",
    ".fa.gz",
    ".cds.fa",
    ".pep.fa",
    ".faa",
    ".fasta",
    ".fna",
    ".fa",
)


def fasta_stem_accession(filename: str) -> str:
    """Accession key from ``<Accession>.faa.gz``-style names (must match BED basename)."""
    fn = filename
    lower = fn.lower()
    for suf in _FASTA_FILENAME_SUFFIXES:
        if lower.endswith(suf):
            return fn[: -len(suf)]
    return os.path.splitext(fn)[0]


def is_sequence_fasta_filename(filename: str) -> bool:
    lower = filename.lower()
    return any(lower.endswith(suf) for suf in _FASTA_FILENAME_SUFFIXES)


@functools.lru_cache(maxsize=1)
def _pigz_executable() -> str | None:
    """Path to ``pigz`` if on PATH; ``None`` means use stdlib :mod:`gzip`."""
    return shutil.which("pigz")


class _PigzTextStream:
    """Context-managed text stream over ``pigz -dc`` (parallel gzip decompress)."""

    __slots__ = ("_proc", "_tw")

    def __init__(self, path: str, *, newline: str | None) -> None:
        pigz = _pigz_executable()
        if not pigz:
            raise RuntimeError("_PigzTextStream requires pigz on PATH")
        self._proc = subprocess.Popen(
            [pigz, "-dc", "--", path],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if self._proc.stdout is None:
            raise OSError("pigz did not provide stdout")
        self._tw = io.TextIOWrapper(
            self._proc.stdout,
            encoding="utf-8",
            errors="replace",
            newline=newline,
        )

    def __enter__(self):
        return self._tw

    def __exit__(self, exc_type, exc, tb):
        try:
            self._tw.close()
        finally:
            err = b""
            if self._proc.stderr:
                err = self._proc.stderr.read() or b""
            rc = self._proc.wait()
            if rc != 0 and err:
                msg = err.decode("utf-8", errors="replace")[:500]
                print(f"pigz -dc failed (exit {rc}): {msg}", file=sys.stderr)
        return False


def _open_gzip_read_text(path: str, *, newline: str | None):
    """Decompress ``.gz`` to text; prefers ``pigz`` when available."""
    if _pigz_executable():
        return _PigzTextStream(path, newline=newline)
    return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline=newline)


def open_text_read(path: str):
    """Text read for plain or gzip (newline='' for csv / line iteration)."""
    if path.lower().endswith(".gz"):
        return _open_gzip_read_text(path, newline="")
    return open(path, "r", newline="", encoding="utf-8", errors="replace", buffering=1 << 20)


def open_fasta_text(path: str):
    """Text read for plain or gzip-compressed FASTA."""
    if path.lower().endswith(".gz"):
        return _open_gzip_read_text(path, newline=None)
    return open(path, "r", newline="", encoding="utf-8", errors="replace", buffering=1 << 20)


class _FastaRec:
    """Minimal duck-type that satisfies ``fasta_header_token_candidates``."""
    __slots__ = ("id", "description", "seq")

    def __init__(self, header: str, seq: str) -> None:
        parts = header.split(None, 1)
        self.id = parts[0] if parts else ""
        self.description = header
        self.seq = seq


def _iter_fasta_fast(handle):
    """Yield ``_FastaRec`` objects without BioPython overhead.

    Roughly 4–6× faster than ``SeqIO.parse`` for simple FASTA files.
    """
    header: str | None = None
    seqbuf: list[str] = []
    for line in handle:
        line = line.rstrip("\n\r")
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield _FastaRec(header, "".join(seqbuf))
            header = line[1:]
            seqbuf = []
        else:
            seqbuf.append(line)
    if header is not None:
        yield _FastaRec(header, "".join(seqbuf))


def _read_fasta_for_db(
    path: str,
    accession: str,
    valid_gene_ids: set,
    is_cds: bool = False,
) -> list[tuple[str, bytes]]:
    """Read one FASTA file in a worker thread; return ``(gene_id, compressed_seq)`` pairs.

    Safe to call from multiple threads since it returns data without touching SQLite.
    Protein sequences have every ``*`` removed (terminal or internal); CDS keeps trailing
    ``*`` strip only, then uppercases.
    """
    results: list[tuple[str, bytes]] = []
    try:
        with open_fasta_text(path) as handle:
            for rec in _iter_fasta_fast(handle):
                gid = resolve_fasta_record_to_db_gene_id(accession, rec, valid_gene_ids)
                if not gid:
                    continue
                if is_cds:
                    seq = rec.seq.upper().rstrip("*")
                else:
                    seq = rec.seq.replace("*", "")
                results.append((gid, zlib.compress(seq.encode("utf-8"), level=6)))
    except (OSError, FileNotFoundError):
        pass
    return results


_BLOB_UNIQ_TABLES = frozenset({"protein_seq_uniq", "cds_seq_uniq"})
_BLOB_MAP_FOR_UNIQ = {"protein_seq_uniq": "protein_seq_map", "cds_seq_uniq": "cds_seq_map"}


def _blob_uniq_id(
    cur: sqlite3.Cursor, cache: dict[bytes, int], uniq_table: str, seq_comp: bytes
) -> int:
    if uniq_table not in _BLOB_UNIQ_TABLES:
        raise ValueError(uniq_table)
    uid = cache.get(seq_comp)
    if uid is not None:
        return uid
    cur.execute(f"INSERT OR IGNORE INTO {uniq_table} (seq_comp) VALUES (?)", (seq_comp,))
    cur.execute(f"SELECT uniq_id FROM {uniq_table} WHERE seq_comp = ?", (seq_comp,))
    row = cur.fetchone()
    if not row:
        raise RuntimeError(f"missing row in {uniq_table} after insert")
    uid = int(row[0])
    cache[seq_comp] = uid
    return uid


def _flush_gene_blob_maps(
    cur: sqlite3.Cursor,
    cache: dict[bytes, int],
    uniq_table: str,
    rows: list[tuple[str, bytes]],
) -> int:
    """Insert (gene_id, uniq_id) into the map table for each (gene_id, seq_comp) row."""
    if not rows:
        return 0
    map_table = _BLOB_MAP_FOR_UNIQ[uniq_table]
    maps: list[tuple[str, int]] = []
    for gid, comp in rows:
        uid = _blob_uniq_id(cur, cache, uniq_table, comp)
        maps.append((gid, uid))
    cur.executemany(
        f"INSERT OR REPLACE INTO {map_table} (gene_id, uniq_id) VALUES (?, ?)",
        maps,
    )
    return len(maps)


def iter_accession_fasta_files(seq_dir: str):
    """Yield ``(path, accession_stem)`` for each sequence file (``.fa``, ``.faa``, ``.gz``, …)."""
    try:
        names = sorted(os.listdir(seq_dir))
    except OSError:
        return
    for fn in names:
        if not is_sequence_fasta_filename(fn):
            continue
        yield os.path.join(seq_dir, fn), fasta_stem_accession(fn)


def normalize_chr_for_coords(chrom: str) -> str:
    """Normalize BED chr strings like 'chr1A' → '1A' for gene_coords."""
    s = str(chrom).strip()
    su = s.upper()
    if su.startswith("CHR"):
        return s[3:]
    return s


# Longest first so ``.bed6.gz`` beats ``.bed.gz`` on odd names, and ``.bed.gz`` beats ``.bed``.
_BED_FILENAME_SUFFIXES = (".bed6.gz", ".bed.gz", ".bed6", ".bed")


def bed_accession_stem(filename: str) -> str | None:
    """Accession key from ``<Accession>.bed.gz``-style names (must match FASTA basename)."""
    lower = filename.lower()
    for suf in _BED_FILENAME_SUFFIXES:
        if lower.endswith(suf):
            return filename[: -len(suf)]
    return None


def load_gene_coords_from_bed(
    cur,
    bed_dir: str,
) -> dict[str, str]:
    """
    Load `input/<dataset>/bed/<Accession>.bed` (or ``.bed.gz``, ``.bed6``, ``.bed6.gz``) into gene_coords.

    BED expected fields (tab-separated):
      chr, start, end, gene_id, <score>, strand, gene_id(again)

    Returns:
      gene_id -> accession mapping (for subsequent pangenome membership loading).
      Gene ids are canonical (no ``accession|`` prefix; see ``canonical_gene_id``).
    """
    inserted = 0
    gene_accession: dict[str, str] = {}

    batch: list[tuple] = []
    batch_size = 20000

    for fn in sorted(os.listdir(bed_dir)):
        accession = bed_accession_stem(fn)
        if not accession:
            continue
        path = os.path.join(bed_dir, fn)
        with open_text_read(path) as f:
            for line in f:
                if not line:
                    continue
                s = line.rstrip("\n")
                if not s or s.startswith("#"):
                    continue
                parts = s.split("\t")
                if len(parts) < 6:
                    continue

                chrom_raw = parts[0]
                chrom = normalize_chr_for_coords(chrom_raw)
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue
                gene_id = canonical_gene_id((parts[3] or "").strip())
                strand = (parts[5] or ".").strip()
                if not gene_id:
                    continue
                subg = subgenome_from_chromosome(chrom)

                gene_accession[gene_id] = accession
                batch.append((gene_id, accession, chrom, start, end, strand, subg))
                inserted += 1

                if len(batch) >= batch_size:
                    cur.executemany(
                        """
                        INSERT OR REPLACE INTO gene_coords
                        (gene_id, accession, chr, start, end, strand, subgenome)
                        VALUES (?,?,?,?,?,?,?)
                        """,
                        batch,
                    )
                    batch = []

    if batch:
        cur.executemany(
            """
            INSERT OR REPLACE INTO gene_coords
            (gene_id, accession, chr, start, end, strand, subgenome)
            VALUES (?,?,?,?,?,?,?)
            """,
            batch,
        )

    print(
        f"  Loaded {inserted:,} BED gene_coord rows into gene_coords (and built accession map)."
    )
    return gene_accession


def iter_pandagma_pan_gene_pairs(pan_membership_tsv_path: str):
    """
    Yield (pan_id, gene_id) from Pandagma pan membership TSV.

    Uses tab-split line parsing (no ``csv`` module) for speed on large clust files.

    Supports:
    - **Wide clust** (e.g. ``18_syn_pan_aug_extra.clust.tsv``): first column is the pan-gene id,
      remaining tab-separated cells are gene IDs (one row per pan-gene cluster).
    - **Hsh** (two columns): ``pan_gene_id<TAB>gene_id`` per row.
    """
    path = os.path.abspath(pan_membership_tsv_path)
    buf = 1 << 20
    with open(path, "r", newline="", buffering=buf) as f:
        first_line = f.readline()
        if not first_line:
            return
        first = first_line.rstrip("\r\n").split("\t")
        wide = len(first) > 2
        if wide:
            pan_id = (first[0] or "").strip()
            for cell in first[1:]:
                gid = canonical_gene_id((cell or "").strip())
                if pan_id and gid:
                    yield (pan_id, gid)
            for line in f:
                if not line or line.lstrip().startswith("#"):
                    continue
                row = line.rstrip("\r\n").split("\t")
                if not row:
                    continue
                pan_id = (row[0] or "").strip()
                for cell in row[1:]:
                    gid = canonical_gene_id((cell or "").strip())
                    if pan_id and gid:
                        yield (pan_id, gid)
        else:
            row = first
            while True:
                if len(row) >= 2:
                    pan_id = (row[0] or "").strip()
                    gid = canonical_gene_id((row[1] or "").strip())
                    if pan_id and gid:
                        yield (pan_id, gid)
                line = f.readline()
                if not line:
                    break
                if line.lstrip().startswith("#"):
                    continue
                row = line.rstrip("\r\n").split("\t")


def _porter6_uniq_id(cur: sqlite3.Cursor, cache: dict[tuple[str, str], int], q3: str, q8: str) -> int:
    key = (q3, q8)
    uid = cache.get(key)
    if uid is not None:
        return uid
    cur.execute("INSERT OR IGNORE INTO porter6_uniq (q3, q8) VALUES (?, ?)", (q3, q8))
    cur.execute(
        "SELECT uniq_id FROM porter6_uniq WHERE q3 = ? AND q8 = ?",
        (q3, q8),
    )
    row = cur.fetchone()
    if not row:
        raise RuntimeError("missing row in porter6_uniq after insert")
    uid = int(row[0])
    cache[key] = uid
    return uid


def load_porter6_from_dir(cur, porter6_dir: str) -> int:
    """
    Load porter6 secondary-structure CSVs into ``porter6_uniq`` + ``porter6_map``.

    Expects paired ``<accession>.q3.csv`` / ``<accession>.q8.csv`` files (``.csv.gz`` ok).
    Columns: ``id`` (gene_id, may carry ``accession|`` prefix), ``q3_decoded`` / ``q8_decoded``.
    Gene IDs are canonicalized via ``canonical_gene_id`` (prefix stripped).
    Rows where *both* q3 and q8 are empty are skipped.
    Returns number of ``porter6_map`` rows written.
    """
    q3_paths: dict[str, str] = {}
    q8_paths: dict[str, str] = {}
    try:
        names = os.listdir(porter6_dir)
    except OSError:
        return 0
    for fn in sorted(names):
        lower = fn.lower()
        if lower.endswith(".q3.csv.gz"):
            stem = fn[: -len(".q3.csv.gz")]
            q3_paths[stem] = os.path.join(porter6_dir, fn)
        elif lower.endswith(".q8.csv.gz"):
            stem = fn[: -len(".q8.csv.gz")]
            q8_paths[stem] = os.path.join(porter6_dir, fn)
        elif lower.endswith(".q3.csv"):
            stem = fn[: -len(".q3.csv")]
            q3_paths[stem] = os.path.join(porter6_dir, fn)
        elif lower.endswith(".q8.csv"):
            stem = fn[: -len(".q8.csv")]
            q8_paths[stem] = os.path.join(porter6_dir, fn)

    def _read_porter6_csv(path: str, col: str) -> dict[str, str]:
        out: dict[str, str] = {}
        try:
            with open_text_read(path) as f:
                for row in csv.DictReader(f):
                    gid = canonical_gene_id((row.get("id") or "").strip())
                    val = (row.get(col) or "").strip()
                    if gid and val:
                        out[gid] = val
        except (OSError, csv.Error):
            pass
        return out

    q3_data: dict[str, str] = {}
    q8_data: dict[str, str] = {}

    all_porter6_tasks = (
        [(path, "q3_decoded") for path in q3_paths.values()] +
        [(path, "q8_decoded") for path in q8_paths.values()]
    )
    if all_porter6_tasks:
        workers = min(8, len(all_porter6_tasks))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futs = {pool.submit(_read_porter6_csv, path, col): col for path, col in all_porter6_tasks}
            for fut in as_completed(futs):
                col = futs[fut]
                result = fut.result()
                if col == "q3_decoded":
                    q3_data.update(result)
                else:
                    q8_data.update(result)

    all_gids = set(q3_data) | set(q8_data)
    batch: list[tuple[str, str, str]] = []
    batch_size = 50_000
    inserted = 0
    uniq_cache: dict[tuple[str, str], int] = {}

    for gid in all_gids:
        q3 = q3_data.get(gid, "")
        q8 = q8_data.get(gid, "")
        if not q3 and not q8:
            continue
        batch.append((gid, q3, q8))
        if len(batch) >= batch_size:
            map_rows: list[tuple[str, int]] = []
            for g, q3v, q8v in batch:
                uid = _porter6_uniq_id(cur, uniq_cache, q3v, q8v)
                map_rows.append((g, uid))
            cur.executemany(
                "INSERT OR REPLACE INTO porter6_map (gene_id, uniq_id) VALUES (?, ?)",
                map_rows,
            )
            inserted += len(map_rows)
            batch = []

    if batch:
        map_rows = []
        for g, q3v, q8v in batch:
            uid = _porter6_uniq_id(cur, uniq_cache, q3v, q8v)
            map_rows.append((g, uid))
        cur.executemany(
            "INSERT OR REPLACE INTO porter6_map (gene_id, uniq_id) VALUES (?, ?)",
            map_rows,
        )
        inserted += len(map_rows)

    return inserted


def build_pandagma_index(
    pan_membership_tsv: str,
    bed_dir: str,
    out_db_path: str,
    *,
    force: bool = False,
    species_id: str = "",
    prot_dir: str | None = None,
    cds_dir: str | None = None,
    porter6_dir: str | None = None,
) -> tuple[bool, dict[str, int] | None]:
    """
    Build a SQLite index for Pandagma pan-gene clusters:
      - ``genes``: member gene id, pan-gene id (``pangene``), accession
      - ``pangene_info``: one row per pan-gene id with member counts
      - ``gene_coords``: from accession-level BED files
      - ``protein_seq_uniq`` / ``protein_seq_map`` and ``cds_seq_uniq`` / ``cds_seq_map``:
        deduplicated zlib-compressed FASTA payloads
      - ``porter6_uniq`` / ``porter6_map``: deduplicated Porter6 q3/q8 strings per gene

    ``species_id`` (folder name under ``input/``, e.g. ``wheat``, ``barley``) controls stored pan-gene
    ids: wheat → ``Traes_<pandagma_pan>``, barley → ``HORVU_<pandagma_pan>``, otherwise unchanged.
    """
    pan_membership_tsv = os.path.abspath(pan_membership_tsv)
    bed_dir = os.path.abspath(bed_dir)
    out_db_path = os.path.abspath(out_db_path)
    prot_dir = os.path.abspath(prot_dir) if prot_dir else ""
    cds_dir = os.path.abspath(cds_dir) if cds_dir else ""

    if not os.path.isfile(pan_membership_tsv):
        raise FileNotFoundError(f"Pandagma pan membership TSV not found: {pan_membership_tsv}")
    if not os.path.isdir(bed_dir):
        raise FileNotFoundError(f"Pandagma bed dir not found: {bed_dir}")

    os.makedirs(os.path.dirname(out_db_path), exist_ok=True)

    if os.path.exists(out_db_path):
        if not force:
            print(f"Database already exists at {out_db_path}")
            print("Run again with --force to overwrite.")
            return False, None
        os.unlink(out_db_path)

    print(f"Building Pandagma pan-gene index from {os.path.basename(pan_membership_tsv)} ...")
    start = time.time()

    conn = sqlite3.connect(out_db_path)
    cur = conn.cursor()
    cur.execute("PRAGMA journal_mode=WAL")
    cur.execute("PRAGMA synchronous=OFF")
    cur.execute("PRAGMA cache_size=-200000")  # ~200 MiB page cache while building
    cur.execute("PRAGMA temp_store=MEMORY")

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS genes (
            gene_id TEXT PRIMARY KEY,
            pangene TEXT NOT NULL,
            accession TEXT NOT NULL
        )
        """
    )

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS pangene_info (
            pangene TEXT PRIMARY KEY,
            clade TEXT,
            gene_count INTEGER
        )
        """
    )

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS gene_coords (
            gene_id TEXT NOT NULL,
            accession TEXT NOT NULL,
            chr TEXT NOT NULL,
            start INTEGER NOT NULL,
            end INTEGER NOT NULL,
            strand TEXT NOT NULL,
            subgenome TEXT,
            chrom_index INTEGER,
            PRIMARY KEY (gene_id, accession)
        )
        """
    )

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS protein_seq_uniq (
            uniq_id INTEGER PRIMARY KEY,
            seq_comp BLOB NOT NULL UNIQUE
        )
        """
    )
    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS protein_seq_map (
            gene_id TEXT PRIMARY KEY,
            uniq_id INTEGER NOT NULL
        )
        """
    )

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS cds_seq_uniq (
            uniq_id INTEGER PRIMARY KEY,
            seq_comp BLOB NOT NULL UNIQUE
        )
        """
    )
    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS cds_seq_map (
            gene_id TEXT PRIMARY KEY,
            uniq_id INTEGER NOT NULL
        )
        """
    )

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS porter6_uniq (
            uniq_id INTEGER PRIMARY KEY,
            q3 TEXT NOT NULL,
            q8 TEXT NOT NULL,
            UNIQUE (q3, q8)
        )
        """
    )
    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS porter6_map (
            gene_id TEXT PRIMARY KEY,
            uniq_id INTEGER NOT NULL
        )
        """
    )

    # 1) Load coordinates + gene→accession mapping.
    if os.path.isdir(bed_dir):
        print(f"  Loading BED gene coords from {bed_dir} ...")
        gid_to_acc = load_gene_coords_from_bed(cur, bed_dir)
    else:
        gid_to_acc = {}

    # 2) Load pan-cluster membership (wide clust.tsv or hsh two-column TSV) into genes.
    gene_batch: list[tuple] = []
    counts: dict[str, int] = defaultdict(int)
    batch_size = 50000
    gene_count = 0

    print("  Loading pan membership from Pandagma TSV ...")
    row_num = 0
    for raw_pan_id, gene_id in iter_pandagma_pan_gene_pairs(pan_membership_tsv):
        row_num += 1
        acc = gid_to_acc.get(gene_id)
        if not acc:
            continue

        pan_id = canonical_pan_id_for_db(raw_pan_id, species_id)
        gene_batch.append((gene_id, pan_id, acc))
        counts[pan_id] += 1
        gene_count += 1

        if len(gene_batch) >= batch_size:
            cur.executemany(
                "INSERT OR IGNORE INTO genes VALUES (?, ?, ?)", gene_batch
            )
            gene_batch = []

        if row_num % 200000 == 0:
            print(f"  Processed {row_num:,} TSV rows, {gene_count:,} genes indexed ...")

    if gene_batch:
        cur.executemany(
            "INSERT OR IGNORE INTO genes VALUES (?, ?, ?)", gene_batch
        )

    # 3) One metadata row per pan-gene cluster id (SQLite table pangene_info).
    pangene_batch = [(hid, "", gcnt) for hid, gcnt in counts.items()]
    cur.executemany(
        "INSERT OR REPLACE INTO pangene_info VALUES (?, ?, ?)", pangene_batch
    )

    # 4) Indexes for speed.
    cur.execute("CREATE INDEX IF NOT EXISTS idx_genes_pangene ON genes(pangene)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_genes_accession ON genes(accession)")

    if counts and os.path.isdir(bed_dir):
        print("  Assigning chrom_index and window index ...")
        assign_chrom_indices(cur)
        cur.execute(
            "CREATE INDEX IF NOT EXISTS idx_coords_accession_chr ON gene_coords(accession, chr, start)"
        )
        cur.execute(
            "CREATE INDEX IF NOT EXISTS idx_coords_gene ON gene_coords(gene_id)"
        )
        cur.execute(
            "CREATE INDEX IF NOT EXISTS idx_coords_chr_window "
            "ON gene_coords(accession, chr, chrom_index)"
        )

    valid_gene_ids = set(r[0] for r in cur.execute("SELECT gene_id FROM genes"))

    # 5) Embed protein sequences (deduplicated blobs + per-gene map).
    # This avoids any FASTA file access / indexing during app runtime.
    protein_uniq_cache: dict[bytes, int] = {}
    if os.path.isdir(prot_dir):
        print(f"  Loading protein sequences from {prot_dir} ...")
        print(f"    Valid gene_id count: {len(valid_gene_ids):,}")
        prot_files = list(iter_accession_fasta_files(prot_dir))
        inserted = 0
        workers = min(8, max(1, len(prot_files)))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futs = {
                pool.submit(_read_fasta_for_db, path, acc, valid_gene_ids, False): path
                for path, acc in prot_files
            }
            batch: list[tuple[str, bytes]] = []
            for fut in as_completed(futs):
                batch.extend(fut.result())
                if len(batch) >= 50_000:
                    inserted += _flush_gene_blob_maps(
                        cur, protein_uniq_cache, "protein_seq_uniq", batch
                    )
                    batch = []
        if batch:
            inserted += _flush_gene_blob_maps(
                cur, protein_uniq_cache, "protein_seq_uniq", batch
            )
        n_prot_uniq = cur.execute("SELECT COUNT(*) FROM protein_seq_uniq").fetchone()[0]
        print(
            f"    Embedded {inserted:,} protein_seq_map rows "
            f"({int(n_prot_uniq):,} distinct protein_seq_uniq)"
        )
    else:
        print(
            f"Warning: prot_dir not found ({prot_dir}); protein tables will be empty.",
            file=sys.stderr,
        )

    cds_uniq_cache: dict[bytes, int] = {}
    if os.path.isdir(cds_dir):
        print(f"  Loading CDS sequences from {cds_dir} ...")
        cds_files = list(iter_accession_fasta_files(cds_dir))
        inserted_cds = 0
        workers_cds = min(8, max(1, len(cds_files)))
        with ThreadPoolExecutor(max_workers=workers_cds) as pool:
            futs_cds = {
                pool.submit(_read_fasta_for_db, path, acc, valid_gene_ids, True): path
                for path, acc in cds_files
            }
            batch_cds: list[tuple[str, bytes]] = []
            for fut in as_completed(futs_cds):
                batch_cds.extend(fut.result())
                if len(batch_cds) >= 50_000:
                    inserted_cds += _flush_gene_blob_maps(
                        cur, cds_uniq_cache, "cds_seq_uniq", batch_cds
                    )
                    batch_cds = []
        if batch_cds:
            inserted_cds += _flush_gene_blob_maps(
                cur, cds_uniq_cache, "cds_seq_uniq", batch_cds
            )
        n_cds_uniq = cur.execute("SELECT COUNT(*) FROM cds_seq_uniq").fetchone()[0]
        print(
            f"    Embedded {inserted_cds:,} cds_seq_map rows "
            f"({int(n_cds_uniq):,} distinct cds_seq_uniq)"
        )
    else:
        print(
            f"Note: cds_dir not found ({cds_dir}); CDS tables left empty (Ka/Ks disabled).",
            file=sys.stderr,
        )

    # 6) Porter6 secondary-structure predictions (optional; skipped if dir absent/empty).
    porter6_dir_abs = os.path.abspath(porter6_dir) if porter6_dir else ""
    if os.path.isdir(porter6_dir_abs):
        print(f"  Loading porter6 secondary structure from {porter6_dir_abs} ...")
        n_p6 = load_porter6_from_dir(cur, porter6_dir_abs)
        if n_p6:
            cur.execute(
                "CREATE INDEX IF NOT EXISTS idx_porter6_map_gene_id ON porter6_map(gene_id)"
            )
        n_p6_uniq = cur.execute("SELECT COUNT(*) FROM porter6_uniq").fetchone()[0]
        print(
            f"    Embedded {n_p6:,} porter6_map rows ({int(n_p6_uniq):,} distinct porter6_uniq)"
        )
    else:
        print(
            f"Note: porter6_dir not found ({porter6_dir_abs}); Porter6 tables left empty.",
            file=sys.stderr,
        )

    stats = stats_from_cursor(cur)
    print(
        f"  Stats: {stats['accessions']} accessions, {stats['pan_genes']:,} pan-genes, "
        f"{stats['genes']:,} genes"
    )
    if stats["genes"] == 0 and gid_to_acc:
        raise RuntimeError(
            "genes table is empty after loading pan membership, but gene_coords has rows. "
            "Pan TSV gene_id values do not match BED gene_coords keys (check ID format, "
            "e.g. HORVU.MOREX.r3 vs HORVU.MOREX.PROJ). Rebuild with a matching pan TSV or "
            "run repair_genes.py --from-coords as a temporary fix."
        )

    conn.commit()
    conn.close()

    elapsed = time.time() - start
    print(f"Done! {len(counts):,} pan-genes indexed in {elapsed:.1f}s")
    print(f"Database saved to {out_db_path}")
    return True, stats


def assign_chrom_indices(cur) -> None:
    """
    For each (accession, chr), set chrom_index = 0..n-1 by genomic order
    (start, end, gene_id). Enables O(window) synteny queries instead of
    loading an entire chromosome per focal gene.
    """
    cur.execute(
        "SELECT DISTINCT accession, chr FROM gene_coords ORDER BY accession, chr"
    )
    pairs = cur.fetchall()
    for row in pairs:
        acc, chrom = row[0], row[1]
        cur.execute(
            """
            SELECT gene_id, accession FROM gene_coords
            WHERE accession = ? AND chr = ?
            ORDER BY start ASC, end ASC, gene_id ASC
            """,
            (acc, chrom),
        )
        genes = cur.fetchall()
        # Use column order (gene_id, accession); works with tuples (default cursor) and sqlite3.Row.
        updates = [(i, g[0], g[1]) for i, g in enumerate(genes)]
        if updates:
            cur.executemany(
                """
                UPDATE gene_coords SET chrom_index = ?
                WHERE gene_id = ? AND accession = ?
                """,
                updates,
            )


def upgrade_gene_coords_chrom_index(db_path: str) -> bool:
    """
    Add chrom_index to an existing DB, populate it, and create the window index.
    Use from Python: ``upgrade_gene_coords_chrom_index("/path/to/species.db")``
    """
    db_path = os.path.abspath(db_path)
    if not os.path.isfile(db_path):
        print(f"Database not found: {db_path}", file=sys.stderr)
        return False
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name='gene_coords'"
    )
    if not cur.fetchone():
        print("No gene_coords table; nothing to do.", file=sys.stderr)
        conn.close()
        return False
    cur.execute("PRAGMA table_info(gene_coords)")
    cols = {r[1] for r in cur.fetchall()}
    if "chrom_index" not in cols:
        cur.execute("ALTER TABLE gene_coords ADD COLUMN chrom_index INTEGER")
        print("Added column gene_coords.chrom_index")
    conn.execute("PRAGMA journal_mode=WAL")
    print("Assigning chrom_index per chromosome (may take a minute) ...")
    t0 = time.time()
    assign_chrom_indices(cur)
    cur.execute(
        "CREATE INDEX IF NOT EXISTS idx_coords_chr_window "
        "ON gene_coords(accession, chr, chrom_index)"
    )
    conn.commit()
    conn.close()
    print(f"Done in {time.time() - t0:.1f}s — synteny queries can use narrow windows.")
    return True


def pick_pandagma_pan_tsv(species_dir: str) -> str | None:
    """Prefer ``*.clust.tsv`` (excluding leftovers / counts); else first ``*.hsh.tsv``."""
    try:
        names = os.listdir(species_dir)
    except OSError:
        return None
    clust = [
        n
        for n in names
        if n.endswith(".clust.tsv")
        and "counts" not in n.lower()
        and "count" not in n.lower()
    ]
    if clust:
        non_lo = [n for n in clust if "leftover" not in n.lower()]
        pick = sorted(non_lo or clust)[0]
        return os.path.join(species_dir, pick)
    hsh = [n for n in names if "hsh" in n.lower() and n.endswith(".tsv")]
    if hsh:
        return os.path.join(species_dir, sorted(hsh)[0])
    return None


def iter_species_input_dirs(input_root: str) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    if not os.path.isdir(input_root):
        return out
    for name in sorted(os.listdir(input_root)):
        if name.startswith("."):
            continue
        path = os.path.join(input_root, name)
        if os.path.isdir(path):
            out.append((name.strip().lower(), path))
    return out


def build_all_from_input_layout(
    repo_root: str,
    *,
    force: bool = False,
    only: list[str] | None = None,
) -> None:
    """
    For each subdirectory of ``input/<species>/``, write ``database/<species>.db``.

    Expects Pandagma/GeneTribe: ``*.clust.tsv`` (or ``*.hsh.tsv``) + ``bed/``
    (``.bed`` / ``.bed6`` / gzip), and optional ``prot/`` and ``cds/`` FASTA directories.
    ``bed`` / ``prot`` may be symlinks to a shared tree (e.g. ``wheat/bed``).

    If ``only`` is set (dataset folder names), build just those and fail loudly if
    any requested name is missing under ``input/``.
    """
    input_root = os.path.join(repo_root, "input")
    db_dir = os.path.join(repo_root, "database")
    os.makedirs(db_dir, exist_ok=True)

    stats_rows: list[tuple[str, dict[str, int]]] = []
    species_dirs = iter_species_input_dirs(input_root)
    if not species_dirs:
        print(f"No species directories under {input_root}", file=sys.stderr)
        return

    if only:
        want = {s.strip().lower() for s in only if s.strip()}
        have = {sid: path for sid, path in species_dirs}
        missing = sorted(want - set(have))
        if missing:
            raise RuntimeError(
                "Requested dataset(s) not found under input/: "
                + ", ".join(missing)
                + f" (have: {', '.join(sorted(have)) or '(none)'})"
            )
        species_dirs = [(sid, have[sid]) for sid in sorted(want)]

    for species_id, sp_dir in species_dirs:
        print(f"\n=== {species_id} ===")
        out_db = os.path.join(db_dir, f"{species_id}.db")
        bed_dir = os.path.join(sp_dir, "bed")
        prot_dir = os.path.join(sp_dir, "prot")
        cds_dir = os.path.join(sp_dir, "cds")
        porter6_dir = os.path.join(sp_dir, "porter6")
        pan_tsv = pick_pandagma_pan_tsv(sp_dir)

        p_prot = prot_dir if os.path.isdir(prot_dir) else None
        p_cds = cds_dir if os.path.isdir(cds_dir) else None
        p_porter6 = porter6_dir if os.path.isdir(porter6_dir) else None

        if pan_tsv and os.path.isdir(bed_dir):
            ok, stats = build_pandagma_index(
                pan_tsv,
                bed_dir,
                out_db,
                force=force,
                species_id=species_id,
                prot_dir=p_prot,
                cds_dir=p_cds,
                porter6_dir=p_porter6,
            )
            if ok and stats:
                stats_rows.append((species_id, stats))
        else:
            print(
                f"  Skip {species_id}: need Pandagma pan TSV (*.clust.tsv or *.hsh.tsv) "
                f"and bed/ under {sp_dir}",
                file=sys.stderr,
            )

    stats_path = os.path.join(db_dir, "stats.tsv")
    if stats_rows:
        merge_write_combined_stats_tsv(stats_path, stats_rows)
        print(f"\nUpdated combined stats: {stats_path}")
    else:
        print("\nNo databases were built; stats.tsv not updated.", file=sys.stderr)


if __name__ == "__main__":
    _root = os.path.dirname(os.path.abspath(__file__))
    _argv = sys.argv[1:]
    if "--chrom-index-only" in _argv:
        print(
            "Chromosome index upgrades are not exposed on the CLI anymore. "
            "Use upgrade_gene_coords_chrom_index(db_path) from a Python shell, "
            "or add a short script if you need this regularly.",
            file=sys.stderr,
        )
        sys.exit(1)
    if "-h" in _argv or "--help" in _argv:
        print(
            "Usage: python build_index.py [--force] [dataset ...]\n"
            "\n"
            "Build database/<dataset>.db from input/<dataset>/ "
            "(*.hsh.tsv or *.clust.tsv + bed/; optional prot/ cds/).\n"
            "With no dataset names, builds every input/*/ folder.\n"
            "Example (GeneTribe wheat): python build_index.py --force wheat_gt",
            file=sys.stderr,
        )
        sys.exit(0)
    _force = "--force" in _argv
    _only = [a for a in _argv if not a.startswith("-")]
    build_all_from_input_layout(
        _root, force=_force, only=_only or None
    )
