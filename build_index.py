"""
Build a SQLite index mapping gene IDs to their HOG and OG groups.

Optionally loads per-accession gene coordinate TSVs into `gene_coords` (wheat GFF-derived).

Optionally loads Porter6 secondary-structure CSVs from `primary/porter6/`: either
`<Accession>.q3.csv` + `<Accession>.q8.csv` (merged by gene id), or legacy `<Accession>.csv`
with id + q3 + q8. Rows are inserted only when `id` matches a gene
already in `genes`. Genes with no Porter row (or accessions not yet loaded) are not stored;
the web app substitutes all-dash (`-`) predictions per residue at display time.
"""
import argparse
import csv
import itertools
import os
import re
import sqlite3
import sys
import time
from collections import defaultdict
import zlib

from Bio import SeqIO

from dataset_stats import stats_from_cursor, write_dataset_stats_tsv


def subgenome_from_chromosome(chrom: str) -> str | None:
    """Infer A/B/D from wheat-style chromosome names like 1A, 4D, 7B."""
    if not chrom:
        return None
    s = str(chrom).strip().upper()
    # BED sometimes uses "chr1A" / "chr4D" naming; coords TSV often uses "1A" / "4D".
    if s.startswith("CHR"):
        s = s[3:]
    # Supported patterns:
    # - "1A" / "4D"  -> A/D at the end
    # - "A1" / "D7"  -> A/D at the start (e.g., chrA1)
    m = re.fullmatch(r"\d+([ABD])", s)
    if m:
        return m.group(1)
    m = re.fullmatch(r"([ABD])\d+", s)
    if m:
        return m.group(1)
    return None


def _porter6_fieldmap(names: list[str] | None) -> dict[str, str]:
    return {k.lower().strip(): k for k in (names or [])}


def _porter6_gene_id(fieldmap: dict[str, str], row: dict[str, str]) -> str:
    for key in ("id", "gene_id"):
        k = fieldmap.get(key)
        if k:
            return (row.get(k) or "").strip()
    return ""


def _porter6_q3_q8_from_row(fieldmap: dict[str, str], row: dict[str, str]) -> tuple[str, str]:
    """q3/q8 from one row; accepts q3 / q3_decoded / q8 / q8_decoded."""
    q3 = ""
    for name in ("q3", "q3_decoded"):
        k = fieldmap.get(name)
        if k:
            q3 = (row.get(k) or "").strip()
            break
    q8 = ""
    for name in ("q8", "q8_decoded"):
        k = fieldmap.get(name)
        if k:
            q8 = (row.get(k) or "").strip()
            break
    return q3, q8


def _porter6_load_split_column(
    path: str, field: str
) -> dict[str, str]:
    """Load id -> sequence string from a single .q3 or .q8 CSV."""
    out: dict[str, str] = {}
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        fieldmap = _porter6_fieldmap(reader.fieldnames)
        if field == "q3":
            col_candidates = ("q3", "q3_decoded")
        else:
            col_candidates = ("q8", "q8_decoded")
        col_key = None
        for name in col_candidates:
            if fieldmap.get(name):
                col_key = fieldmap[name]
                break
        if not col_key:
            return out
        for row in reader:
            gid = _porter6_gene_id(fieldmap, row)
            if not gid:
                continue
            out[gid] = (row.get(col_key) or "").strip()
    return out


def _porter6_flush_batch(
    cur,
    batch: list[tuple[str, str, str, str]],
) -> None:
    if not batch:
        return
    cur.executemany(
        """
        INSERT OR REPLACE INTO porter6_ss
        (gene_id, accession, q3, q8)
        VALUES (?,?,?,?)
        """,
        batch,
    )


def load_porter6_csvs(cur, porter6_dir: str, valid_gene_ids: set[str]) -> int:
    """
    Load Porter6 predictions from ``primary/porter6/``:

    - **Split (preferred):** ``<Accession>.q3.csv`` + ``<Accession>.q8.csv``, merged on gene id.
      Headers may use ``q3`` / ``q3_decoded`` and ``q8`` / ``q8_decoded``.
    - **Legacy:** ``<Accession>.csv`` with id + q3 + q8 (or *_decoded) on one row.

    Only rows whose id is in ``valid_gene_ids`` are inserted (same policy as protein_seqs).

    Returns:
        Number of rows inserted/replaced.
    """
    inserted = 0
    batch: list[tuple[str, str, str, str]] = []
    batch_size = 20000

    def flush() -> None:
        nonlocal batch
        if batch:
            _porter6_flush_batch(cur, batch)
            batch = []

    all_csv = [fn for fn in os.listdir(porter6_dir) if fn.endswith(".csv")]
    bases_q3 = {fn[: -len(".q3.csv")] for fn in all_csv if fn.endswith(".q3.csv")}
    bases_q8 = {fn[: -len(".q8.csv")] for fn in all_csv if fn.endswith(".q8.csv")}
    split_bases = bases_q3 & bases_q8
    unpaired_split = (bases_q3 | bases_q8) - split_bases
    if unpaired_split:
        print(
            f"  Porter6: skip accessions missing .q3 or .q8 pair: {', '.join(sorted(unpaired_split))}",
            file=sys.stderr,
        )

    for accession in sorted(split_bases):
        p3 = os.path.join(porter6_dir, f"{accession}.q3.csv")
        p8 = os.path.join(porter6_dir, f"{accession}.q8.csv")
        q3_map = _porter6_load_split_column(p3, "q3")
        q8_map = _porter6_load_split_column(p8, "q8")
        if not q3_map and not q8_map:
            continue
        ids = set(q3_map) & set(q8_map)
        only_q3 = set(q3_map) - set(q8_map)
        only_q8 = set(q8_map) - set(q3_map)
        if only_q3 or only_q8:
            print(
                f"  Porter6 {accession}: {len(only_q3)} ids only in .q3, "
                f"{len(only_q8)} only in .q8 (using intersection {len(ids):,}).",
                file=sys.stderr,
            )
        for gid in ids:
            if gid not in valid_gene_ids:
                continue
            batch.append((gid, accession, q3_map[gid], q8_map[gid]))
            inserted += 1
            if len(batch) >= batch_size:
                flush()

    flush()

    # Legacy single-file <Accession>.csv (not .q3.csv / .q8.csv)
    legacy_names = [
        fn
        for fn in sorted(all_csv)
        if not fn.endswith(".q3.csv") and not fn.endswith(".q8.csv")
    ]
    for fn in legacy_names:
        accession = os.path.splitext(fn)[0]
        if accession in split_bases:
            continue
        path = os.path.join(porter6_dir, fn)
        with open(path, "r", newline="") as f:
            reader = csv.DictReader(f)
            names = reader.fieldnames or []
            if not names:
                continue
            fieldmap = _porter6_fieldmap(names)
            for row in reader:
                gid = _porter6_gene_id(fieldmap, row)
                if not gid or gid not in valid_gene_ids:
                    continue
                q3, q8 = _porter6_q3_q8_from_row(fieldmap, row)
                batch.append((gid, accession, q3, q8))
                inserted += 1
                if len(batch) >= batch_size:
                    flush()

    flush()
    return inserted


def normalize_chr_for_coords(chrom: str) -> str:
    """Normalize BED chr strings like 'chr1A' → '1A' for gene_coords."""
    s = str(chrom).strip()
    su = s.upper()
    if su.startswith("CHR"):
        return s[3:]
    return s


def load_gene_coords_from_bed(
    cur,
    bed_dir: str,
) -> dict[str, str]:
    """
    Load `input/<dataset>/bed/<Accession>.bed` into gene_coords.

    BED expected fields (tab-separated):
      chr, start, end, gene_id, <score>, strand, gene_id(again)

    Returns:
      gene_id -> accession mapping (for subsequent pangenome membership loading).
    """
    inserted = 0
    gene_accession: dict[str, str] = {}

    batch: list[tuple] = []
    batch_size = 20000

    for fn in sorted(os.listdir(bed_dir)):
        if not fn.endswith(".bed"):
            continue
        path = os.path.join(bed_dir, fn)
        accession = os.path.splitext(fn)[0]
        with open(path, "r", newline="") as f:
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
                gene_id = (parts[3] or "").strip()
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

    Supports:
    - **Wide clust** (e.g. ``18_syn_pan_aug_extra.clust.tsv``): first column is pan_id,
      remaining tab-separated cells are gene IDs (one row per pan).
    - **Hsh** (two columns): ``pan_id<TAB>gene_id`` per row.
    """
    path = os.path.abspath(pan_membership_tsv_path)
    with open(path, "r", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        first = next(reader, None)
        if not first:
            return
        wide = len(first) > 2
        if wide:
            pan_id = (first[0] or "").strip()
            for cell in first[1:]:
                gid = (cell or "").strip()
                if pan_id and gid:
                    yield (pan_id, gid)
            for row in reader:
                if not row:
                    continue
                pan_id = (row[0] or "").strip()
                for cell in row[1:]:
                    gid = (cell or "").strip()
                    if pan_id and gid:
                        yield (pan_id, gid)
        else:
            for row in itertools.chain([first], reader):
                if len(row) < 2:
                    continue
                pan_id = (row[0] or "").strip()
                gid = (row[1] or "").strip()
                if pan_id and gid:
                    yield (pan_id, gid)


def build_pandagma_index(
    pan_membership_tsv: str,
    bed_dir: str,
    out_db_path: str,
    *,
    force: bool = False,
    og_constant: str = "pan",
    prot_dir: str | None = None,
    cds_dir: str | None = None,
    porter6_dir: str | None = None,
) -> bool:
    """
    Build a SQLite index for Pandagma pan-clusters:
      - genes table: gene_id → pan_id, og(=pan), accession
      - hog_info table: pan_id → pan info (og/clade placeholders)
      - gene_coords table: from accession-level BED files
      - protein_seqs table: gene_id → compressed protein sequence (no FASTA reads later)
      - porter6_ss table: gene_id, accession → q3/q8 secondary structure strings
    """
    pan_membership_tsv = os.path.abspath(pan_membership_tsv)
    bed_dir = os.path.abspath(bed_dir)
    out_db_path = os.path.abspath(out_db_path)
    if prot_dir is None:
        prot_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "primary", "prot")
    prot_dir = os.path.abspath(prot_dir)
    if cds_dir is None:
        cds_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "primary", "cds")
    cds_dir = os.path.abspath(cds_dir) if cds_dir else ""
    if porter6_dir is None:
        porter6_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "primary", "porter6")
    porter6_dir = os.path.abspath(porter6_dir) if porter6_dir else ""

    if not os.path.isfile(pan_membership_tsv):
        raise FileNotFoundError(f"Pandagma pan membership TSV not found: {pan_membership_tsv}")
    if not os.path.isdir(bed_dir):
        raise FileNotFoundError(f"Pandagma bed dir not found: {bed_dir}")

    os.makedirs(os.path.dirname(out_db_path), exist_ok=True)

    if os.path.exists(out_db_path):
        if not force:
            print(f"Database already exists at {out_db_path}")
            print("Run again with --force to overwrite.")
            return False
        os.unlink(out_db_path)

    print(f"Building Pandagma pan index from {os.path.basename(pan_membership_tsv)} ...")
    start = time.time()

    conn = sqlite3.connect(out_db_path)
    cur = conn.cursor()
    cur.execute("PRAGMA journal_mode=WAL")
    cur.execute("PRAGMA synchronous=OFF")

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS genes (
            gene_id TEXT PRIMARY KEY,
            hog TEXT NOT NULL,
            og TEXT NOT NULL,
            accession TEXT NOT NULL
        )
        """
    )

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS hog_info (
            hog TEXT PRIMARY KEY,
            og TEXT NOT NULL,
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
        CREATE TABLE IF NOT EXISTS protein_seqs (
            gene_id TEXT PRIMARY KEY,
            seq_comp BLOB NOT NULL
        )
        """
    )

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS cds_seqs (
            gene_id TEXT PRIMARY KEY,
            seq_comp BLOB NOT NULL
        )
        """
    )

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS porter6_ss (
            gene_id TEXT NOT NULL,
            accession TEXT NOT NULL,
            q3 TEXT NOT NULL,
            q8 TEXT NOT NULL,
            PRIMARY KEY (gene_id, accession)
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
    for pan_id, gene_id in iter_pandagma_pan_gene_pairs(pan_membership_tsv):
        row_num += 1
        acc = gid_to_acc.get(gene_id)
        if not acc:
            continue

        gene_batch.append((gene_id, pan_id, og_constant, acc))
        counts[pan_id] += 1
        gene_count += 1

        if len(gene_batch) >= batch_size:
            cur.executemany(
                "INSERT OR IGNORE INTO genes VALUES (?, ?, ?, ?)", gene_batch
            )
            gene_batch = []

        if row_num % 200000 == 0:
            print(f"  Processed {row_num:,} pan–gene pairs, {gene_count:,} genes ...")

    if gene_batch:
        cur.executemany(
            "INSERT OR IGNORE INTO genes VALUES (?, ?, ?, ?)", gene_batch
        )

    # 3) hog_info (pan_id metadata).
    hog_batch = [(hid, og_constant, "", gcnt) for hid, gcnt in counts.items()]
    cur.executemany(
        "INSERT OR REPLACE INTO hog_info VALUES (?, ?, ?, ?)", hog_batch
    )

    # 4) Indexes for speed.
    cur.execute("CREATE INDEX IF NOT EXISTS idx_genes_hog ON genes(hog)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_genes_og ON genes(og)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_genes_accession ON genes(accession)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_hog_og ON hog_info(og)")

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

    # 5) Embed protein sequences directly into SQLite.
    # This avoids any FASTA file access / indexing during app runtime.
    if os.path.isdir(prot_dir):
        print(f"  Loading protein sequences from {prot_dir} ...")
        valid_gene_ids = set(r[0] for r in cur.execute("SELECT gene_id FROM genes"))
        print(f"    Valid gene_id count: {len(valid_gene_ids):,}")

        batch: list[tuple[str, bytes]] = []
        batch_size = 20000
        inserted = 0

        def flush():
            nonlocal batch, inserted
            if not batch:
                return
            cur.executemany(
                "INSERT OR REPLACE INTO protein_seqs (gene_id, seq_comp) VALUES (?, ?)",
                batch,
            )
            inserted += len(batch)
            batch = []

        # Parse each accession-level FASTA once and insert only sequences present in genes.
        for fn in sorted(os.listdir(prot_dir)):
            if not (fn.endswith(".fa") or fn.endswith(".fasta") or fn.endswith(".fna")):
                continue
            fasta_path = os.path.join(prot_dir, fn)
            try:
                for record in SeqIO.parse(fasta_path, "fasta"):
                    gid = (record.id or "").strip()
                    if gid not in valid_gene_ids:
                        continue
                    clean_seq = str(record.seq).rstrip("*")
                    comp = zlib.compress(clean_seq.encode("utf-8"), level=6)
                    batch.append((gid, comp))
                    if len(batch) >= batch_size:
                        flush()
            except FileNotFoundError:
                continue
            # Flush remaining records between files so memory stays bounded.
            flush()

        flush()
        print(f"    Embedded {inserted:,} protein sequences into protein_seqs")
    else:
        print(f"Warning: prot_dir not found ({prot_dir}); protein_seqs will be empty.", file=sys.stderr)

    if os.path.isdir(cds_dir):
        print(f"  Loading CDS sequences from {cds_dir} ...")
        batch_cds: list[tuple[str, bytes]] = []
        inserted_cds = 0
        batch_size_cds = 20000

        def flush_cds():
            nonlocal batch_cds, inserted_cds
            if not batch_cds:
                return
            cur.executemany(
                "INSERT OR REPLACE INTO cds_seqs (gene_id, seq_comp) VALUES (?, ?)",
                batch_cds,
            )
            inserted_cds += len(batch_cds)
            batch_cds = []

        for fn in sorted(os.listdir(cds_dir)):
            if not (fn.endswith(".fa") or fn.endswith(".fasta") or fn.endswith(".fna")):
                continue
            fasta_path = os.path.join(cds_dir, fn)
            try:
                for record in SeqIO.parse(fasta_path, "fasta"):
                    gid = (record.id or "").strip()
                    if gid not in valid_gene_ids:
                        continue
                    clean_seq = str(record.seq).upper().rstrip("*")
                    comp = zlib.compress(clean_seq.encode("utf-8"), level=6)
                    batch_cds.append((gid, comp))
                    if len(batch_cds) >= batch_size_cds:
                        flush_cds()
            except FileNotFoundError:
                continue
            flush_cds()

        flush_cds()
        print(f"    Embedded {inserted_cds:,} CDS sequences into cds_seqs")
    else:
        print(
            f"Note: cds_dir not found ({cds_dir}); cds_seqs left empty (Ka/Ks disabled).",
            file=sys.stderr,
        )

    if os.path.isdir(porter6_dir):
        print(f"  Loading Porter6 secondary structure from {porter6_dir} ...")
        t_p6 = time.time()
        gid_set = set(r[0] for r in cur.execute("SELECT gene_id FROM genes"))
        n_p6 = load_porter6_csvs(cur, porter6_dir, gid_set)
        print(f"    Loaded {n_p6:,} porter6_ss rows in {time.time() - t_p6:.1f}s")
        # gene_id-only index for lookups without accession (PK is gene_id+accession)
        cur.execute(
            "CREATE INDEX IF NOT EXISTS idx_porter6_ss_gene_id ON porter6_ss(gene_id)"
        )
        cur.execute(
            "CREATE INDEX IF NOT EXISTS idx_porter6_ss_accession ON porter6_ss(accession)"
        )
    else:
        print(
            f"Note: porter6_dir not found ({porter6_dir}); porter6_ss left empty.",
            file=sys.stderr,
        )

    stats = stats_from_cursor(cur)
    stats_tsv = os.path.join(
        os.path.dirname(os.path.abspath(out_db_path)), "dataset_stats.tsv"
    )
    write_dataset_stats_tsv(stats, stats_tsv)
    print(
        f"Wrote {stats_tsv} ({stats['accessions']} acc, {stats['ogs']} OGs, {stats['hogs']:,} HOGs/pans, {stats['genes']:,} genes)"
    )

    conn.commit()
    conn.close()

    elapsed = time.time() - start
    print(f"Done! {len(counts):,} pans indexed in {elapsed:.1f}s")
    print(f"Database saved to {out_db_path}")
    return True


def load_gene_coords(cur, coords_dir: str) -> int:
    """
    Load `input/<species>/coords/<Accession>.tsv` files.

    Expected header: gene_id, chr, start, end, strand
    """
    inserted = 0
    batch: list[tuple] = []
    batch_size = 20000

    for fn in sorted(os.listdir(coords_dir)):
        if not fn.endswith(".tsv"):
            continue
        path = os.path.join(coords_dir, fn)
        accession = os.path.splitext(fn)[0]
        with open(path, "r", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            # normalize keys
            fieldmap = {k.lower().strip(): k for k in (reader.fieldnames or [])}

            def col(row, name):
                k = fieldmap.get(name.lower())
                return (row.get(k) or "").strip() if k else ""

            for row in reader:
                gid = col(row, "gene_id")
                chrom = col(row, "chr")
                if not gid or not chrom:
                    continue
                try:
                    start = int(col(row, "start"))
                    end = int(col(row, "end"))
                except ValueError:
                    continue
                strand = col(row, "strand") or "."
                sg = subgenome_from_chromosome(chrom)
                batch.append((gid, accession, chrom, start, end, strand, sg))
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
            batch = []

    return inserted


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
    Use: python build_index.py --chrom-index-only --out input/wheat/panwheat.db
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


def build_index(
    n0_path: str,
    out_db_path: str,
    *,
    force: bool = False,
    coords_dir: str | None = None,
    porter6_dir: str | None = None,
) -> bool:
    """
    Build one SQLite index. Returns True if a new DB was written, False if skipped
    (output exists and force=False).
    """
    n0_path = os.path.abspath(n0_path)
    out_db_path = os.path.abspath(out_db_path)
    if porter6_dir is None:
        porter6_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "primary", "porter6")
    porter6_dir = os.path.abspath(porter6_dir) if porter6_dir else ""

    if not os.path.exists(n0_path):
        raise FileNotFoundError(f"N0 TSV not found: {n0_path}")

    os.makedirs(os.path.dirname(out_db_path), exist_ok=True)

    if os.path.exists(out_db_path):
        if not force:
            print(f"Database already exists at {out_db_path}")
            print("Run again with --force to overwrite.")
            return False
        os.unlink(out_db_path)

    print(f"Building gene index from {os.path.basename(n0_path)} ...")
    start = time.time()

    conn = sqlite3.connect(out_db_path)
    cur = conn.cursor()

    cur.execute("PRAGMA journal_mode=WAL")
    cur.execute("PRAGMA synchronous=OFF")

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS genes (
            gene_id TEXT PRIMARY KEY,
            hog TEXT NOT NULL,
            og TEXT NOT NULL,
            accession TEXT NOT NULL
        )
        """
    )

    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS hog_info (
            hog TEXT PRIMARY KEY,
            og TEXT NOT NULL,
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
        CREATE TABLE IF NOT EXISTS porter6_ss (
            gene_id TEXT NOT NULL,
            accession TEXT NOT NULL,
            q3 TEXT NOT NULL,
            q8 TEXT NOT NULL,
            PRIMARY KEY (gene_id, accession)
        )
        """
    )

    gene_count = 0
    hog_count = 0

    with open(n0_path, "r") as f:
        try:
            csv.field_size_limit(sys.maxsize)
        except OverflowError:
            csv.field_size_limit(2147483647)

        reader = csv.reader(f, delimiter="\t")
        header = next(reader)

        accession_cols = []
        for i, col in enumerate(header[3:], start=3):
            acc_name = col.replace(".primary.protein", "")
            accession_cols.append((i, acc_name))

        batch = []
        hog_batch = []

        for row_num, row in enumerate(reader, start=2):
            if len(row) < 3:
                continue
            hog = row[0]
            og = row[1]
            clade = row[2]
            row_gene_count = 0

            for col_idx, acc_name in accession_cols:
                if col_idx >= len(row):
                    continue
                cell = row[col_idx].strip()
                if not cell:
                    continue
                for gene_id in cell.split(", "):
                    gene_id = gene_id.strip()
                    if gene_id:
                        batch.append((gene_id, hog, og, acc_name))
                        gene_count += 1
                        row_gene_count += 1

            hog_batch.append((hog, og, clade, row_gene_count))
            hog_count += 1

            if len(batch) >= 50000:
                cur.executemany(
                    "INSERT OR IGNORE INTO genes VALUES (?, ?, ?, ?)", batch
                )
                batch = []

            if len(hog_batch) >= 10000:
                cur.executemany(
                    "INSERT OR IGNORE INTO hog_info VALUES (?, ?, ?, ?)", hog_batch
                )
                hog_batch = []

            if row_num % 20000 == 0:
                print(f"  Processed {row_num:,} rows, {gene_count:,} genes ...")

    if batch:
        cur.executemany("INSERT OR IGNORE INTO genes VALUES (?, ?, ?, ?)", batch)
    if hog_batch:
        cur.executemany(
            "INSERT OR IGNORE INTO hog_info VALUES (?, ?, ?, ?)", hog_batch
        )

    print("Creating indexes ...")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_genes_hog ON genes(hog)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_genes_og ON genes(og)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_genes_accession ON genes(accession)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_hog_og ON hog_info(og)")

    if coords_dir:
        coords_dir = os.path.abspath(coords_dir)
        if os.path.isdir(coords_dir):
            print(f"Loading coordinates from {coords_dir} ...")
            t0 = time.time()
            ncoord = load_gene_coords(cur, coords_dir)
            print(
                f"  Inserted/updated {ncoord:,} coordinate rows in {time.time() - t0:.1f}s"
            )
            print("  Assigning chrom_index for fast synteny windows ...")
            t1 = time.time()
            assign_chrom_indices(cur)
            print(f"    chrom_index done in {time.time() - t1:.1f}s")
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
        else:
            print(f"Warning: --coords-dir not found: {coords_dir}", file=sys.stderr)

    if os.path.isdir(porter6_dir):
        print(f"Loading Porter6 secondary structure from {porter6_dir} ...")
        t_p6 = time.time()
        gid_set = set(r[0] for r in cur.execute("SELECT gene_id FROM genes"))
        n_p6 = load_porter6_csvs(cur, porter6_dir, gid_set)
        print(f"  Inserted {n_p6:,} porter6_ss rows in {time.time() - t_p6:.1f}s")
        # gene_id-only index for lookups without accession (PK is gene_id+accession)
        cur.execute(
            "CREATE INDEX IF NOT EXISTS idx_porter6_ss_gene_id ON porter6_ss(gene_id)"
        )
        cur.execute(
            "CREATE INDEX IF NOT EXISTS idx_porter6_ss_accession ON porter6_ss(accession)"
        )
    else:
        print(
            f"Note: porter6_dir not found ({porter6_dir}); porter6_ss left empty.",
            file=sys.stderr,
        )

    stats = stats_from_cursor(cur)
    stats_tsv = os.path.join(os.path.dirname(os.path.abspath(out_db_path)), "dataset_stats.tsv")
    write_dataset_stats_tsv(stats, stats_tsv)
    print(f"Wrote {stats_tsv} ({stats['accessions']} acc, {stats['ogs']} OGs, {stats['hogs']:,} HOGs, {stats['genes']:,} genes)")

    conn.commit()
    conn.close()

    elapsed = time.time() - start
    print(f"Done! {hog_count:,} HOGs, {gene_count:,} genes indexed in {elapsed:.1f}s")
    print(f"Database saved to {out_db_path}")
    return True


def run_default_barley_then_wheat(
    root: str, *, force: bool, porter6_dir: str | None = None
) -> None:
    """Barley from OrthoFinder N0 (optional); wheat from Pandagma clust + BED (preferred) or N0 fallback."""
    barley_n0 = os.path.join(root, "input", "barley", "BPGv2_N0.tsv")
    barley_out = os.path.join(root, "input", "barley", "panbarley.db")
    wheat_n0 = os.path.join(root, "input", "wheat", "N0.tsv")
    wheat_out = os.path.join(root, "input", "wheat", "panwheat.db")
    wheat_coords = os.path.join(root, "input", "wheat", "coords")
    wheat_clust = os.path.join(root, "primary", "18_syn_pan_aug_extra.clust.tsv")
    wheat_bed = os.path.join(root, "primary", "bed")
    wheat_pandagma_out = os.path.join(root, "input", "wheat", "panwheat_pandagma.db")

    print("=== Barley (no coordinates) ===")
    if not os.path.isfile(barley_n0):
        print(f"  Skipping barley: N0 not found: {barley_n0}")
    else:
        build_index(barley_n0, barley_out, force=force, coords_dir=None, porter6_dir=porter6_dir)

    print("\n=== Wheat (Pandagma: clust.tsv + primary/bed) ===")
    if os.path.isfile(wheat_clust) and os.path.isdir(wheat_bed):
        build_pandagma_index(
            wheat_clust,
            wheat_bed,
            wheat_pandagma_out,
            force=force,
            porter6_dir=porter6_dir,
        )
    elif os.path.isfile(wheat_n0):
        print(
            f"  Note: Pandagma inputs not found ({wheat_clust} and/or {wheat_bed}); "
            "falling back to OrthoFinder N0 → panwheat.db.",
            file=sys.stderr,
        )
        coords_arg = wheat_coords if os.path.isdir(wheat_coords) else None
        if coords_arg is None:
            print(
                f"  Warning: coords directory not found ({wheat_coords}); "
                "building wheat DB without gene_coords.",
                file=sys.stderr,
            )
        build_index(wheat_n0, wheat_out, force=force, coords_dir=coords_arg, porter6_dir=porter6_dir)
    else:
        print(
            f"  Skipping wheat: need either Pandagma ({wheat_clust} + {wheat_bed}) "
            f"or OrthoFinder N0 ({wheat_n0}).",
            file=sys.stderr,
        )


def resolve_single_build_paths(args, root: str) -> tuple[str, str, str | None]:
    """Defaults when user passes any of --n0 / --out / --coords-dir."""
    default_barley_n0 = os.path.join(root, "input", "barley", "BPGv2_N0.tsv")
    default_barley_out = os.path.join(root, "input", "barley", "panbarley.db")

    n0 = args.n0_path or default_barley_n0
    if args.out_db_path:
        out = args.out_db_path
    elif args.n0_path:
        d = os.path.dirname(os.path.abspath(args.n0_path))
        dnorm = d.replace("\\", "/").lower()
        if "wheat" in dnorm:
            out = os.path.join(d, "panwheat.db")
        else:
            out = os.path.join(d, "panbarley.db")
    else:
        out = default_barley_out
    coords = args.coords_dir
    return n0, out, coords


if __name__ == "__main__":
    _root = os.path.dirname(os.path.abspath(__file__))

    parser = argparse.ArgumentParser(
        description=(
            "Build SQLite pan-genome indexes: barley from OrthoFinder N0; wheat from "
            "Pandagma primary/18_syn_pan_aug_extra.clust.tsv + primary/bed (default), "
            "or OrthoFinder N0 if Pandagma inputs are missing."
        ),
        epilog=(
            "With no --n0, --out, or --coords-dir: builds barley (N0) then wheat "
            "(Pandagma clust + primary/bed → input/wheat/panwheat_pandagma.db).\n"
            "Examples:\n"
            "  Both:    python build_index.py\n"
            "  Rebuild: python build_index.py --force\n"
            "  N0-only: python build_index.py --n0 input/wheat/N0.tsv "
            "--out input/wheat/panwheat.db --coords-dir input/wheat/coords --force\n"
            "  Pandagma CLI: python build_index.py --pandagma-pan-tsv primary/18_syn_pan_aug_extra.clust.tsv "
            "--pandagma-bed-dir primary/bed --out input/wheat/panwheat_pandagma.db --force"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--n0",
        dest="n0_path",
        default=None,
        help="Input N0 TSV (HOG/OG membership table).",
    )
    parser.add_argument(
        "--out",
        dest="out_db_path",
        default=None,
        help="Output SQLite database path.",
    )
    parser.add_argument(
        "--coords-dir",
        dest="coords_dir",
        default=None,
        help="Optional directory of <Accession>.tsv files with gene_id,chr,start,end,strand.",
    )
    parser.add_argument(
        "--pandagma-pan-tsv",
        "--pandagma-hsh",
        dest="pandagma_pan_tsv",
        default=None,
        help=(
            "Pandagma pan membership: wide clust.tsv (pan_id + gene columns) or "
            "two-column hsh (pan_id<TAB>gene_id per row)."
        ),
    )
    parser.add_argument(
        "--pandagma-bed-dir",
        dest="pandagma_bed_dir",
        default=None,
        help="Directory of accession BED files used to populate gene_coords.",
    )
    parser.add_argument(
        "--pandagma-og-constant",
        dest="pandagma_og_constant",
        default="pan",
        help="Constant OG label to store in DB (Pandagma output has no OG IDs).",
    )
    parser.add_argument(
        "--pandagma-prot-dir",
        dest="pandagma_prot_dir",
        default=None,
        help="Directory of accession protein FASTA files (default: <repo>/primary/prot).",
    )
    parser.add_argument(
        "--pandagma-cds-dir",
        dest="pandagma_cds_dir",
        default=None,
        help="Directory of accession CDS FASTA files for Ka/Ks (default: <repo>/primary/cds).",
    )
    parser.add_argument(
        "--porter6-dir",
        dest="porter6_dir",
        default=None,
        help=(
            "Directory of <Accession>.csv Porter6 outputs (columns id,q3,q8). "
            "Default: <repo>/primary/porter6 if present."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite existing output DB if it exists.",
    )
    parser.add_argument(
        "--chrom-index-only",
        action="store_true",
        help=(
            "Only add/populate gene_coords.chrom_index and the window index "
            "(use on an existing wheat DB; requires --out)."
        ),
    )
    args = parser.parse_args()

    if args.pandagma_pan_tsv and not args.pandagma_bed_dir:
        _default_bed = os.path.join(_root, "primary", "bed")
        if os.path.isdir(_default_bed):
            args.pandagma_bed_dir = _default_bed

    if args.chrom_index_only:
        if not args.out_db_path:
            print("--chrom-index-only requires --out <path/to/panwheat.db>", file=sys.stderr)
            sys.exit(1)
        ok = upgrade_gene_coords_chrom_index(args.out_db_path)
        sys.exit(0 if ok else 1)

    # Pandagma mode: build pan index directly from clust/hsh TSV + BED.
    if args.pandagma_pan_tsv and args.pandagma_bed_dir:
        out_db = (
            args.out_db_path
            if args.out_db_path
            else os.path.join(_root, "input", "wheat", "panwheat_pandagma.db")
        )
        ok = build_pandagma_index(
            args.pandagma_pan_tsv,
            args.pandagma_bed_dir,
            out_db,
            force=args.force,
            og_constant=args.pandagma_og_constant,
            prot_dir=args.pandagma_prot_dir,
            cds_dir=args.pandagma_cds_dir,
            porter6_dir=args.porter6_dir,
        )
        sys.exit(0 if ok else 1)
    elif args.pandagma_pan_tsv or args.pandagma_bed_dir:
        print(
            "Pandagma mode requires both --pandagma-pan-tsv and --pandagma-bed-dir.",
            file=sys.stderr,
        )
        sys.exit(1)

    custom = (
        args.n0_path is not None
        or args.out_db_path is not None
        or args.coords_dir is not None
    )

    if not custom:
        run_default_barley_then_wheat(_root, force=args.force, porter6_dir=args.porter6_dir)
    else:
        n0, out, coords = resolve_single_build_paths(args, _root)
        build_index(
            n0,
            out,
            force=args.force,
            coords_dir=coords,
            porter6_dir=args.porter6_dir,
        )
