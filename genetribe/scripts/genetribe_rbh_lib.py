#!/usr/bin/env python3
"""Shared helpers for GeneTribe RBH → pan-gene clustering."""

from __future__ import annotations

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
                    raise RuntimeError(
                        f"bad BED line (need >=4 cols) in {bed}: {line[:80]!r}"
                    )
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
        f"  Finalize only after all pair jobs succeed. Sample:\n{hint}"
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
