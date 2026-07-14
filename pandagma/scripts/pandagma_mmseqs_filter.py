#!/usr/bin/env python3
"""
Expand MMseqs representative/member rows to allowed member/member pairs.

Input cluster file format (Pandagma MMseqs output):
  chrX__GeneA__start__end__strand<TAB>chrY__GeneB__start__end__strand

This script groups rows by representative (column 1), then emits all-vs-all
pairs inside each representative's member set, filtered by allowed chromosome
pairs from the ``expected_chr_matches`` include. Ingest BED molecules must be
``chr`` + the token from that file, case-insensitive: token ``1A`` → ``chr1A``,
``2H`` → ``chr2H``. The script strips a leading ``chr`` and compares the rest to
the token (case-insensitive), then checks the pair against the include.

Output format matches Pandagma filter output (8 tab-separated columns):
  chr1 gene1 start1 end1 chr2 gene2 start2 end2

Use a genome directory as the sole argument (see below) to process every
coordinate-style ``*_cluster.tsv`` under ``<genome>/work/03_mmseqs/`` and write
``<genome>/work/04_dag/<stem>_matches.tsv`` for each.

Important for ``pandagma pan -s dagchainer``: each work file must be named exactly
``04_dag/<Genome1>.x.<Genome2>_matches.tsv`` with a single ``.x.`` in the stem.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path


def parse_expected_chr_matches(inc_path: Path) -> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()
    for raw in inc_path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        if s in {"expected_chr_matches=(", ")"}:
            continue
        fields = s.split()
        if len(fields) >= 2:
            a, b = fields[0], fields[1]
            pairs.add((a, b))
            pairs.add((b, a))
    return pairs


def _token_spelling_by_lower(allowed: set[tuple[str, str]]) -> dict[str, str]:
    """Map ``token.lower()`` -> spelling as in the include (for ``(t1,t2) in allowed``)."""
    out: dict[str, str] = {}
    for a, b in allowed:
        for t in (a, b):
            out.setdefault(t.lower(), t)
    return out


def chrom_to_config_token(chrom: str, token_by_lower: dict[str, str]) -> str | None:
    """
    BED molecule is ``f'chr{{token}}'`` (any case). Strip ``chr``, match remainder to token.
    """
    s = chrom.strip().lower()
    if len(s) >= 3 and s.startswith("chr"):
        s = s[3:]
    return token_by_lower.get(s)


def split_hashed(g: str) -> tuple[str, str, str, str]:
    """
    Return (chrom, gene_id, start, end) from:
      chr4B__Traes...__486939488__486944776__+
    or without strand field.
    """
    parts = g.strip().split("__")
    if len(parts) < 4:
        raise ValueError(f"Unexpected hashed gene field: {g}")
    chrom, gene_id, start, end = parts[0], parts[1], parts[2], parts[3]
    return chrom, gene_id, start, end


def edge_key_eight_cols(cols: list[str]) -> tuple[str, ...]:
    if len(cols) < 8:
        raise ValueError("Need 8 columns")
    a, b = cols[:4], cols[4:8]
    return tuple(a + b) if a <= b else tuple(b + a)


def is_coordinate_cluster_format(cluster_path: Path, sample_lines: int = 50) -> bool:
    """True if file looks like Pandagma positional cluster TSV (hashed headers with __)."""
    n = 0
    for raw in cluster_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not raw.strip():
            continue
        cols = raw.split("\t")
        if len(cols) < 2:
            continue
        if "__" not in cols[0] or "__" not in cols[1]:
            return False
        n += 1
        if n >= sample_lines:
            break
    return n > 0


def expand_cluster_members_to_matches(
    cluster: Path,
    out: Path,
    allowed: set[tuple[str, str]],
    merge_existing: bool,
) -> int:
    """Build deduplicated 8-column match rows; return number of rows written."""
    by_rep: dict[str, list[str]] = defaultdict(list)

    for line in cluster.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip():
            continue
        cols = line.split("\t")
        if len(cols) < 2:
            continue
        rep, member = cols[0], cols[1]
        by_rep[rep].append(rep)
        by_rep[rep].append(member)

    out.parent.mkdir(parents=True, exist_ok=True)
    by_key: dict[tuple[str, ...], str] = {}

    if merge_existing and out.is_file():
        for line in out.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            by_key[edge_key_eight_cols(cols)] = line.rstrip("\n")

    token_by_lower = _token_spelling_by_lower(allowed)
    for rep, genes in by_rep.items():
        uniq = sorted(set(genes))
        for a, b in combinations(uniq, 2):
            c1, g1, s1, e1 = split_hashed(a)
            c2, g2, s2, e2 = split_hashed(b)
            t1 = chrom_to_config_token(c1, token_by_lower)
            t2 = chrom_to_config_token(c2, token_by_lower)
            if t1 is None or t2 is None or (t1, t2) not in allowed:
                continue
            line8 = "\t".join([c1, g1, s1, e1, c2, g2, s2, e2])
            cols = line8.split("\t")
            by_key[edge_key_eight_cols(cols)] = line8

    out.write_text("\n".join(by_key.values()) + ("\n" if by_key else ""), encoding="utf-8")
    return len(by_key)


def find_expected_inc(genome_dir: Path) -> Path:
    cfg = genome_dir / "config"
    for pat in ("*_expected_chr_matches.inc", "*expected_chr_matches*.inc"):
        found = sorted(cfg.glob(pat))
        if found:
            return found[0]
    raise FileNotFoundError(f"No expected_chr_matches include under {cfg}")


def run_all_comparisons(
    cluster_anchor: Path,
    allowed: set[tuple[str, str]],
    merge_existing: bool,
) -> int:
    mmseqs_dir = cluster_anchor.parent
    dag_dir = mmseqs_dir.parent / "04_dag"
    cluster_files = sorted(mmseqs_dir.glob("*_cluster.tsv"))
    if not cluster_files:
        print(f"No *_cluster.tsv files in {mmseqs_dir}", flush=True)
        return 1
    bad_neighbors = sorted(dag_dir.glob("*_members_matches.tsv")) if dag_dir.is_dir() else []
    if bad_neighbors:
        print(
            "WARNING: remove these before dagchainer (wrong stem parses as extra genome): "
            + ", ".join(str(p) for p in bad_neighbors),
            flush=True,
        )
    grand_total = 0
    for cluster in cluster_files:
        if not is_coordinate_cluster_format(cluster):
            print(f"Skipping (not coordinate-hashed format): {cluster.name}", flush=True)
            continue
        stem = (
            cluster.name[: -len("_cluster.tsv")]
            if cluster.name.endswith("_cluster.tsv")
            else cluster.stem
        )
        out = dag_dir / f"{stem}_matches.tsv"
        n = expand_cluster_members_to_matches(cluster, out, allowed, merge_existing)
        grand_total += n
        print(f"Wrote {n} deduplicated 8-column rows to: {out}", flush=True)
    print(f"Total rows (all comparisons): {grand_total}", flush=True)
    return 0


def main() -> int:
    repo_root = Path(__file__).resolve().parent.parent
    ap = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    ap.add_argument(
        "genome",
        nargs="?",
        default=None,
        help=(
            "Genome/analysis directory (absolute path or name under repo root). "
            "Uses <genome>/work/03_mmseqs and <genome>/config/*expected_chr_matches*.inc"
        ),
    )
    ap.add_argument(
        "--cluster-tsv",
        type=Path,
        default=None,
        help="MMseqs representative/member TSV (implies directory for multi-file mode when used with --all-comparisons)",
    )
    ap.add_argument(
        "--expected-inc",
        type=Path,
        default=None,
        help="Path to expected_chr_matches include file",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output 8-column matches (single-file mode only; default: work/04_dag/<stem>_matches.tsv)",
    )
    ap.add_argument(
        "--merge-existing",
        action="store_true",
        help="If output exists, keep its 8-column rows and add expanded pairs (deduplicated)",
    )
    ap.add_argument(
        "--all-comparisons",
        action="store_true",
        help="Process every *_cluster.tsv in the same directory as --cluster-tsv",
    )
    ap.add_argument(
        "--no-all-comparisons",
        action="store_true",
        help="With genome argument, only process --cluster-tsv instead of all cluster files",
    )
    args = ap.parse_args()

    genome_dir: Path | None = None
    if args.genome:
        gd = Path(args.genome)
        if not gd.is_absolute():
            gd = repo_root / gd
        genome_dir = gd.resolve()
        if not genome_dir.is_dir():
            print(f"ERROR: genome directory not found: {genome_dir}", file=sys.stderr, flush=True)
            return 1
        mmseqs_dir = genome_dir / "work" / "03_mmseqs"
        if args.expected_inc is None:
            try:
                args.expected_inc = find_expected_inc(genome_dir)
            except FileNotFoundError as e:
                print(f"ERROR: {e}", file=sys.stderr, flush=True)
                return 1
        if args.cluster_tsv is None:
            clusters = sorted(mmseqs_dir.glob("*_cluster.tsv"))
            if not clusters:
                print(f"ERROR: no *_cluster.tsv in {mmseqs_dir}", file=sys.stderr, flush=True)
                return 1
            args.cluster_tsv = clusters[0]
        if not args.no_all_comparisons:
            args.all_comparisons = True

    if args.expected_inc is None:
        print("ERROR: set --expected-inc or pass a genome directory.", file=sys.stderr, flush=True)
        return 1
    if args.cluster_tsv is None:
        print("ERROR: set --cluster-tsv or pass a genome directory.", file=sys.stderr, flush=True)
        return 1

    cluster_anchor = args.cluster_tsv.resolve()
    allowed = parse_expected_chr_matches(args.expected_inc.resolve())

    if args.all_comparisons:
        return run_all_comparisons(cluster_anchor, allowed, args.merge_existing)

    stem = cluster_anchor.name[: -len("_cluster.tsv")] if cluster_anchor.name.endswith("_cluster.tsv") else cluster_anchor.stem
    out = args.out.resolve() if args.out else cluster_anchor.parent.parent / "04_dag" / f"{stem}_matches.tsv"

    bad_neighbors = sorted(out.parent.glob("*_members_matches.tsv"))
    if bad_neighbors:
        print(
            "WARNING: remove these before dagchainer (wrong stem parses as extra genome): "
            + ", ".join(str(p) for p in bad_neighbors),
            flush=True,
        )

    n = expand_cluster_members_to_matches(cluster_anchor, out, allowed, args.merge_existing)
    print(f"Wrote {n} deduplicated 8-column rows to: {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
