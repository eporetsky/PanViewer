#!/usr/bin/env python3
"""
Build pan-genes from GeneTribe pairwise ``*.RBH`` files.

Production rule (subgenome-aware + collinear):
  1. Keep an RBH only when both genes share a chromosome *group* across
     homeologous subgenomes (wheat ``1A+1B+1D``, oat ``1A+1C+1D``,
     barley ``1H``…``7H``).
  2. Keep that edge only if the gene pair also sits in a collinear block
     (``*.collinearity_info`` / ``*.colinearity_info``, jcvi ``*.anchors``,
     or ``*.block_pos``).

Writes under ``<analysis>/work/``:
  - genetribe_subgenome_pans.hsh.tsv    pan_id<TAB>gene_id
  - genetribe_subgenome_pans.clust.tsv  wide cluster rows

Genes never seen in a kept RBH become size-1 pans (from BED gene IDs).
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from genetribe_rbh_lib import (  # noqa: E402
    UnionFind,
    add_bed_singletons,
    assign_pan_ids,
    count_accessions,
    find_rbh_files,
    print_component_report,
)

_CHROM_RE = re.compile(r"^(?:chr)?(\d+)([ABCDH])$", re.I)


@dataclass(frozen=True)
class GeneLocus:
    accession: str
    chrom: str
    chrom_group: str  # "1"…"7" — homeologous group across subgenomes
    index: int
    start: int
    end: int


def chrom_group_from_chrom(chrom: str) -> str | None:
    m = _CHROM_RE.fullmatch((chrom or "").strip())
    return m.group(1) if m else None


def load_gene_loci(analysis: Path) -> dict[str, GeneLocus]:
    acc_root = analysis / "accessions"
    if not acc_root.is_dir():
        raise RuntimeError(f"missing {acc_root}")
    by_chrom: dict[tuple[str, str], list[tuple[int, int, int, str, str]]] = defaultdict(
        list
    )
    n_bed = 0
    n_skip = 0
    for acc_dir in sorted(p for p in acc_root.iterdir() if p.is_dir()):
        bed = acc_dir / f"{acc_dir.name}.bed"
        if not bed.is_file():
            continue
        n_bed += 1
        acc = acc_dir.name
        with bed.open(encoding="utf-8", errors="replace") as fh:
            for line_i, line in enumerate(fh):
                if not line.strip() or line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) < 4:
                    continue
                chrom, gid = cols[0].strip(), cols[3].strip()
                if not gid:
                    continue
                group = chrom_group_from_chrom(chrom)
                if group is None:
                    n_skip += 1
                    continue
                try:
                    start = int(cols[1])
                except ValueError:
                    start = line_i
                try:
                    end = int(cols[2])
                except (ValueError, IndexError):
                    end = start
                exact = f"{group}{chrom.strip()[-1].upper()}"
                by_chrom[(acc, exact)].append((start, end, line_i, gid, group))
    loci: dict[str, GeneLocus] = {}
    counts: dict[str, int] = defaultdict(int)
    for (acc, exact), rows in by_chrom.items():
        rows.sort(key=lambda t: (t[0], t[1], t[3]))
        for idx, (start, end, _li, gid, group) in enumerate(rows):
            prev = loci.get(gid)
            if prev is not None and (
                prev.chrom_group != group or prev.accession != acc
            ):
                raise RuntimeError(
                    f"gene {gid} mapped to both {prev.accession}:{prev.chrom} "
                    f"and {acc}:{exact}"
                )
            loci[gid] = GeneLocus(
                accession=acc,
                chrom=exact,
                chrom_group=group,
                index=idx,
                start=start,
                end=end,
            )
            counts[group] += 1
    if not loci:
        raise RuntimeError(
            f"no genes with N+[ABCDH] chromosomes under {acc_root} "
            "(expected names like 1A, 2B, 7H, chr1A, …)"
        )
    print(
        f"subgenome chrom groups: {len(loci):,} genes from {n_bed} accessions "
        f"(group counts={dict(sorted(counts.items(), key=lambda kv: int(kv[0])))}; "
        f"unparsed chrom: {n_skip:,})",
        file=sys.stderr,
    )
    return loci


def chrom_key(chrom: str) -> str | None:
    group = chrom_group_from_chrom(chrom)
    if group is None:
        return None
    return f"{group}{chrom.strip()[-1].upper()}"


def _pair_key(a: str, b: str) -> tuple[str, str] | None:
    a, b = a.strip(), b.strip()
    if not a or not b or a == b:
        return None
    return (a, b) if a < b else (b, a)


def iter_colinearity_pairs(path: Path):
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            for pair in cols[5].split(";"):
                pair = pair.strip()
                if "," not in pair:
                    continue
                a, b = pair.split(",", 1)
                yield a, b


def first_collinearity_file(pair_dir: Path) -> Path | None:
    files = sorted(pair_dir.glob("*.collinearity_info")) + sorted(
        pair_dir.glob("*.colinearity_info")
    )
    return files[0] if files else None


def parse_block_pos(path: Path) -> list[tuple[str, int, int, str, int, int]]:
    blocks: list[tuple[str, int, int, str, int, int]] = []
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split("\t") if "\t" in line else line.split()
            if len(cols) < 6:
                continue
            c1, c2 = chrom_key(cols[0]), chrom_key(cols[3])
            if c1 is None or c2 is None:
                continue
            try:
                s1, e1, s2, e2 = int(cols[1]), int(cols[2]), int(cols[4]), int(cols[5])
            except ValueError:
                continue
            if e1 < s1:
                s1, e1 = e1, s1
            if e2 < s2:
                s2, e2 = e2, s2
            blocks.append((c1, s1, e1, c2, s2, e2))
            if c1 != c2 or s1 != s2 or e1 != e2:
                blocks.append((c2, s2, e2, c1, s1, e1))
    return blocks


def _in_interval(pos: int, lo: int, hi: int) -> bool:
    return lo <= pos <= hi


def rbh_in_collinear_block(
    la: GeneLocus, lb: GeneLocus, blocks: list[tuple[str, int, int, str, int, int]]
) -> bool:
    mid_a = (la.start + la.end) // 2
    mid_b = (lb.start + lb.end) // 2
    for c1, s1, e1, c2, s2, e2 in blocks:
        if (
            la.chrom == c1
            and lb.chrom == c2
            and _in_interval(mid_a, s1, e1)
            and _in_interval(mid_b, s2, e2)
        ):
            return True
    return False


def collinear_keys_for_rbh(
    pair_dir: Path, rbh_keys: set[tuple[str, str]]
) -> tuple[str, set[tuple[str, str]] | None, list]:
    info = first_collinearity_file(pair_dir)
    if info is not None:
        hits: set[tuple[str, str]] = set()
        for a, b in iter_colinearity_pairs(info):
            key = _pair_key(a, b)
            if key is not None and key in rbh_keys:
                hits.add(key)
        return info.name, hits, []
    go = pair_dir / "genetribe_output"
    hits = set()
    src = ""
    if go.is_dir():
        for path in sorted(go.glob("*.anchors"))[:1]:
            src = f"genetribe_output/{path.name}"
            with path.open(encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    s = line.strip()
                    if not s or s.startswith("#") or s.startswith("###"):
                        continue
                    cols = s.split("\t") if "\t" in s else s.split()
                    if len(cols) < 2:
                        continue
                    key = _pair_key(cols[0], cols[1])
                    if key is not None and key in rbh_keys:
                        hits.add(key)
    if hits:
        return src, hits, []
    blocks: list[tuple[str, int, int, str, int, int]] = []
    sources: list[str] = []
    for path in sorted(pair_dir.glob("*.block_pos")):
        got = parse_block_pos(path)
        if got:
            blocks.extend(got)
            sources.append(path.name)
    if blocks:
        return ",".join(sources), None, blocks
    raise RuntimeError(
        f"no collinearity files in {pair_dir} "
        "(need *.collinearity_info, genetribe_output/*.anchors, or *.block_pos)"
    )


def _pct(n: int, d: int) -> str:
    return f"{(100.0 * n / d):.1f}%" if d else "n/a"


_WORKER_LOCI: dict[str, GeneLocus] = {}


def _scan_one_rbh_file(path_str: str) -> dict:
    loci = _WORKER_LOCI
    path = Path(path_str)
    rows: list[tuple[str, str, GeneLocus, GeneLocus]] = []
    n_edges = n_missing = n_cross = n_same_exact = n_homoeolog = 0
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            a, b = cols[0].strip(), cols[1].strip()
            if not a or not b or a == b:
                continue
            n_edges += 1
            la, lb = loci.get(a), loci.get(b)
            if la is None or lb is None:
                n_missing += 1
                continue
            if la.chrom_group != lb.chrom_group:
                n_cross += 1
                continue
            if la.chrom == lb.chrom:
                n_same_exact += 1
            else:
                n_homoeolog += 1
            rows.append((a, b, la, lb))

    rbh_keys = {_pair_key(a, b) for a, b, _la, _lb in rows}
    rbh_keys.discard(None)
    src, col_keys, col_blocks = collinear_keys_for_rbh(path.parent, rbh_keys)
    kept: list[tuple[str, str]] = []
    n_noncol = 0
    for a, b, la, lb in rows:
        key = _pair_key(a, b)
        ok = False
        if col_keys is not None:
            ok = key in col_keys
        elif col_blocks:
            ok = rbh_in_collinear_block(la, lb, col_blocks)
        if ok:
            kept.append((a, b))
        else:
            n_noncol += 1
    return {
        "n_edges": n_edges,
        "n_missing": n_missing,
        "n_cross": n_cross,
        "n_same_exact": n_same_exact,
        "n_homoeolog": n_homoeolog,
        "n_noncol": n_noncol,
        "col_src": src,
        "kept": kept,
    }


def cluster_subgenome_collinear(
    pairs_root: Path,
    loci: dict[str, GeneLocus],
    *,
    workers: int | None = None,
) -> dict[str, list[str]]:
    global _WORKER_LOCI
    rbh_files = find_rbh_files(pairs_root)
    uf = UnionFind()
    if workers is None:
        workers = int(os.environ.get("SLURM_CPUS_PER_TASK") or "1")
    workers = max(1, workers)
    _WORKER_LOCI = loci

    n_edges = n_kept = n_cross = n_missing = 0
    n_same_exact = n_homoeolog = n_noncol = 0
    col_sources: dict[str, int] = defaultdict(int)
    n_files = len(rbh_files)
    done = 0

    def _apply(result: dict) -> None:
        nonlocal n_edges, n_kept, n_cross, n_missing
        nonlocal n_same_exact, n_homoeolog, n_noncol, done
        n_edges += result["n_edges"]
        n_missing += result["n_missing"]
        n_cross += result["n_cross"]
        n_same_exact += result["n_same_exact"]
        n_homoeolog += result["n_homoeolog"]
        n_noncol += result["n_noncol"]
        if result["col_src"]:
            col_sources[result["col_src"]] += 1
        for a, b in result["kept"]:
            uf.union(a, b)
            n_kept += 1
        done += 1
        if done == 1 or done % 50 == 0 or done == n_files:
            print(f"  scanned {done}/{n_files} RBH files", file=sys.stderr, flush=True)

    print(
        f"scanning {n_files} RBH files with {workers} worker(s) "
        f"(subgenome-aware + collinear)",
        file=sys.stderr,
        flush=True,
    )
    if workers == 1 or n_files < 4:
        for path in rbh_files:
            _apply(_scan_one_rbh_file(str(path)))
    else:
        try:
            from multiprocessing import get_context

            ctx = get_context("fork")
            pool_kw = {"max_workers": workers, "mp_context": ctx}
        except (ValueError, ImportError):
            pool_kw = {"max_workers": workers}
        with ProcessPoolExecutor(**pool_kw) as pool:
            futs = [pool.submit(_scan_one_rbh_file, str(p)) for p in rbh_files]
            for fut in as_completed(futs):
                _apply(fut.result())

    n_same_group = n_same_exact + n_homoeolog
    print(
        f"RBH files: {len(rbh_files)}; edges: {n_edges:,}; "
        f"missing coords: {n_missing:,}; "
        f"dropped different chrom group: {n_cross:,} ({_pct(n_cross, n_edges)}); "
        f"same chrom group: {n_same_group:,} "
        f"(exact {n_same_exact:,}; homeologue {n_homoeolog:,}); "
        f"dropped non-collinear: {n_noncol:,} ({_pct(n_noncol, n_same_group)}); "
        f"kept: {n_kept:,}",
        file=sys.stderr,
    )
    return uf.components()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("analysis", type=Path, help="Directory with accessions/ and work/")
    ap.add_argument("--work", type=Path, default=None)
    ap.add_argument("--prefix", default="OGI", help="Pan ID prefix (default OGI)")
    ap.add_argument(
        "--no-singletons",
        action="store_true",
        help="Do not add size-1 pans for genes absent from kept RBHs",
    )
    args = ap.parse_args()

    analysis = args.analysis.resolve()
    work = (args.work or (analysis / "work")).resolve()
    pairs_root = work / "pairs"
    n_acc = count_accessions(analysis)

    try:
        loci = load_gene_loci(analysis)
        components = cluster_subgenome_collinear(pairs_root, loci)
        if not args.no_singletons:
            components = add_bed_singletons(components, analysis)
    except RuntimeError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    print(f"accessions={n_acc}", file=sys.stderr)
    print_component_report(components)
    pans = assign_pan_ids(components, args.prefix.strip("_") or "OGI")

    work.mkdir(parents=True, exist_ok=True)
    hsh_path = work / "genetribe_subgenome_pans.hsh.tsv"
    clust_path = work / "genetribe_subgenome_pans.clust.tsv"
    with hsh_path.open("w", encoding="utf-8") as hsh, clust_path.open(
        "w", encoding="utf-8"
    ) as clust:
        for pan_id, members in pans:
            clust.write(pan_id + "\t" + "\t".join(members) + "\n")
            for gid in members:
                hsh.write(f"{pan_id}\t{gid}\n")

    n_genes = sum(len(m) for _, m in pans)
    print(f"Wrote {len(pans)} pans, {n_genes:,} gene memberships")
    print(f"Wrote {hsh_path}")
    print(f"Wrote {clust_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
