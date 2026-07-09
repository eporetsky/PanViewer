"""Fast protein track ordering for Collinearity (k-mer Jaccard + UPGMA).

Designed for ~50–200 aligned sequences: O(n² × L) to hash + O(n³) UPGMA with small n.
No scipy/sourmash subprocess; pure Python for predictable latency.
"""

from __future__ import annotations

import zlib


def _kmer_hashes(seq: str, k: int, mod: int) -> set[int]:
    s = "".join(c for c in seq.upper() if c.isalpha())
    if len(s) < k:
        return set()
    out: set[int] = set()
    for i in range(len(s) - k + 1):
        km = s[i : i + k]
        x = zlib.adler32(km.encode("ascii", "ignore")) & 0xFFFFFFFF
        out.add(x % mod)
    return out


def _jaccard(a: set[int], b: set[int]) -> float:
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    inter = len(a & b)
    union = len(a | b)
    return inter / union if union else 0.0


def _upgma_leaf_order(dist: list[list[float]], n0: int) -> list[int]:
    """UPGMA (average linkage) on distance matrix; return leaf indices left-to-right."""
    if n0 <= 1:
        return list(range(n0))

    nmax = 2 * n0
    d = [[1e18] * nmax for _ in range(nmax)]
    for i in range(n0):
        for j in range(n0):
            d[i][j] = dist[i][j]

    sz = [1] * n0 + [0] * (nmax - n0)
    child: list[tuple[int, int] | None] = [None] * nmax
    nxt = n0
    alive = list(range(n0))

    while len(alive) > 1:
        bi = bj = -1
        bd = 1e18
        for ii in range(len(alive)):
            i = alive[ii]
            for j in alive[ii + 1 :]:
                if d[i][j] < bd:
                    bd, bi, bj = d[i][j], i, j
        assert bi >= 0 and bj >= 0
        k = nxt
        nxt += 1
        si, sj = sz[bi], sz[bj]
        sk = si + sj
        sz[k] = sk
        child[k] = (bi, bj)
        for m in alive:
            if m == bi or m == bj:
                continue
            val = (si * d[bi][m] + sj * d[bj][m]) / sk
            d[k][m] = d[m][k] = val
        alive = [m for m in alive if m not in (bi, bj)] + [k]

    root = alive[0]

    def leaves(u: int) -> list[int]:
        if u < n0:
            return [u]
        pair = child[u]
        if pair is None:
            return [u]
        a, b = pair
        return leaves(a) + leaves(b)

    return leaves(root)


def protein_track_order_permutation(
    gene_ids: list[str],
    sequences: dict[str, str],
    *,
    k: int = 5,
    mod: int = 4096,
) -> list[int]:
    """Return a permutation of range(len(gene_ids)) — UPGMA leaf order under 1 - Jaccard."""
    n = len(gene_ids)
    if n <= 1:
        return list(range(n))
    bags: list[set[int]] = []
    for gid in gene_ids:
        seq = (sequences.get(gid) or "").strip()
        bags.append(_kmer_hashes(seq, k=k, mod=mod) if seq else set())

    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            jac = _jaccard(bags[i], bags[j])
            d = 1.0 - jac
            dist[i][j] = dist[j][i] = d

    order = _upgma_leaf_order(dist, n)
    if len(order) != n or set(order) != set(range(n)):
        return list(range(n))
    return order
