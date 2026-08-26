#!/usr/bin/env python3
"""
Build PanViewer pan-gene membership from GeneTribe pairwise ``*.RBH`` files.

This is the production clustering used for wheat/oat/barley ``*.genetribe.db``:
connected components on Reciprocal Best Hit edges (RGI Ortholog Gene Index).

GeneTribe itself does **not** emit pan-genes. Pair jobs write ``.RBH`` / ``.SBH``
/ ``.one2many``; this script only reads ``*.RBH``.

Writes (same names as the existing Atlas finalize):
  - genetribe_pans.hsh.tsv   pan_id<TAB>gene_id
  - genetribe_pans.clust.tsv wide clust rows

Genes absent from any RBH edge become size-1 pans (BED gene IDs) unless
``--no-singletons``.

A separate MCL experiment lives in ``genetribe_mcl_to_pans.py`` and writes
``genetribe_mcl_pans.*`` so it cannot overwrite these files.
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}
        self.rank: dict[str, int] = {}

    def add(self, x: str) -> None:
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0

    def find(self, x: str) -> str:
        self.add(x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1

    def components(self) -> dict[str, list[str]]:
        groups: dict[str, list[str]] = defaultdict(list)
        for x in self.parent:
            groups[self.find(x)].append(x)
        return groups


def parse_rbh_file(path: Path) -> list[tuple[str, str]]:
    """Parse GeneTribe *.RBH (geneA, geneB, optional extra columns)."""
    edges: list[tuple[str, str]] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                cols = line.split()
            if len(cols) < 2:
                continue
            a, b = cols[0].strip(), cols[1].strip()
            if a and b:
                edges.append((a, b))
    return edges


def load_genes_from_beds(analysis: Path) -> set[str]:
    genes: set[str] = set()
    acc_root = analysis / "accessions"
    if not acc_root.is_dir():
        raise RuntimeError(f"missing accessions dir: {acc_root}")
    for bed in sorted(acc_root.glob("*/*.bed")):
        with bed.open("r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) < 4:
                    raise RuntimeError(f"bad BED line (need >=4 cols) in {bed}: {line[:80]!r}")
                gid = cols[3].strip()
                if gid:
                    genes.add(gid)
    return genes


def find_rbh_files(pairs_root: Path) -> list[Path]:
    if not pairs_root.is_dir():
        raise RuntimeError(f"missing pairs directory: {pairs_root}")
    files = sorted(pairs_root.glob("*/*.RBH"))
    if files:
        return files
    pair_dirs = sorted(p for p in pairs_root.iterdir() if p.is_dir())
    n_done = sum(1 for p in pair_dirs if (p / ".done").is_file())
    n_log = sum(1 for p in pair_dirs if (p / "pair.log").is_file())
    sample = []
    for p in pair_dirs[:5]:
        names = sorted(x.name for x in p.iterdir() if x.is_file())[:12]
        sample.append(f"  {p.name}: files={names or '(empty)'}")
    hint = "\n".join(sample) if sample else "  (no pair subdirectories)"
    raise RuntimeError(
        f"no *.RBH files under {pairs_root}\n"
        f"  pair dirs: {len(pair_dirs)}; .done markers: {n_done}; pair.log: {n_log}\n"
        f"  Finalize only after all pair jobs succeed. Check:\n"
        f"    ls {pairs_root} | head\n"
        f"    find {pairs_root} -name '*.RBH' | head\n"
        f"  Sample pair dirs:\n{hint}"
    )


def count_accessions(analysis: Path) -> int:
    acc_root = analysis / "accessions"
    if not acc_root.is_dir():
        return 0
    n = 0
    for d in acc_root.iterdir():
        if d.is_dir() and (d / f"{d.name}.bed").is_file():
            n += 1
    return n


SIZE_BUCKETS = (
    (1, 1, "1"),
    (2, 10, "2-10"),
    (11, 50, "11-50"),
    (51, 200, "51-200"),
    (201, 500, "201-500"),
    (501, 5000, "501-5000"),
    (5001, 10**18, "5000+"),
)


def print_component_report(components: dict[str, list[str]]) -> None:
    sized = sorted(
        ((len(m), sorted(m)) for m in components.values()),
        key=lambda t: (-t[0], t[1][0] if t[1] else ""),
    )
    n_genes = sum(n for n, _ in sized)
    print(f"components: {len(sized)}; genes in graph: {n_genes}")
    print("size histogram (n_pans, n_genes):")
    for lo, hi, label in SIZE_BUCKETS:
        pans = [(n, m) for n, m in sized if lo <= n <= hi]
        if not pans:
            continue
        print(f"  {label:>8}: {len(pans):>8} pans  {sum(n for n, _ in pans):>10,} genes")
    print("largest pans:")
    for n, members in sized[:12]:
        sample = ", ".join(members[:3])
        extra = "…" if n > 3 else ""
        print(f"  {n:>8} genes  e.g. {sample}{extra}")


def assign_pan_ids(
    components: dict[str, list[str]],
    prefix: str,
) -> list[tuple[str, list[str]]]:
    """
    Stable pan IDs: sort each component's genes, sort components by first gene,
    then number OGI_000001 ...
    """
    comps = [sorted(members) for members in components.values()]
    comps.sort(key=lambda members: (members[0], len(members), members))
    out: list[tuple[str, list[str]]] = []
    for i, members in enumerate(comps, start=1):
        pan_id = f"{prefix}_{i:06d}"
        out.append((pan_id, members))
    return out


def add_bed_singletons(
    components: dict[str, list[str]], analysis: Path
) -> dict[str, list[str]]:
    uf = UnionFind()
    for members in components.values():
        if not members:
            continue
        first = members[0]
        uf.add(first)
        for g in members[1:]:
            uf.union(first, g)
    all_genes = load_genes_from_beds(analysis)
    in_graph = set(uf.parent)
    extra = 0
    for g in all_genes:
        if g not in in_graph:
            uf.add(g)
            extra += 1
    print(f"singletons from BED (not in graph): {extra:,}")
    return uf.components()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "analysis",
        type=Path,
        help="Analysis directory (contains accessions/ and work/)",
    )
    ap.add_argument(
        "--work",
        type=Path,
        default=None,
        help="Work directory (default: <analysis>/work)",
    )
    ap.add_argument(
        "--prefix",
        default="OGI",
        help="Pan-gene ID prefix (default: OGI → OGI_000001)",
    )
    ap.add_argument(
        "--no-singletons",
        action="store_true",
        help="Do not add size-1 pans for genes absent from the RBH graph",
    )
    args = ap.parse_args()

    analysis = args.analysis.resolve()
    work = (args.work or (analysis / "work")).resolve()
    pairs_root = work / "pairs"
    n_acc = count_accessions(analysis)

    try:
        rbh_files = find_rbh_files(pairs_root)
        uf = UnionFind()
        n_edges = 0
        for path in rbh_files:
            for a, b in parse_rbh_file(path):
                uf.union(a, b)
                n_edges += 1
        print(f"RBH files: {len(rbh_files)}; edges: {n_edges}; accessions: {n_acc}")
        components = uf.components()
        if not args.no_singletons:
            components = add_bed_singletons(components, analysis)
    except RuntimeError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1

    print_component_report(components)
    pans = assign_pan_ids(components, args.prefix.strip("_") or "OGI")

    work.mkdir(parents=True, exist_ok=True)
    hsh_path = work / "genetribe_pans.hsh.tsv"
    clust_path = work / "genetribe_pans.clust.tsv"

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
