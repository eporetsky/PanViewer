"""
Single normalization for gene identifiers across BED, pan TSV, FASTA, and the app.

Pandagma-style ``{accession}|{locus}`` headers and BED name fields are reduced to
plain ``{locus}`` so SQLite ``genes.gene_id`` / ``gene_coords.gene_id`` / sequence
tables never mix prefixed and unprefixed forms.
"""


def fasta_header_token_candidates(record) -> list[str]:
    """
    Ordered distinct tokens from a Biopython ``SeqRecord`` header (id + description words).

    Headers like ``>MorexV3|HORVU.MOREX... description`` are covered by ``id``; if the
    pipe-form lives only after the first space, later tokens are tried too.
    """
    rid = str(getattr(record, "id", "") or "").strip()
    desc = str(getattr(record, "description", "") or "").strip()
    parts = desc.split() if desc else []
    out: list[str] = []
    seen: set[str] = set()
    for t in [rid] + [p for p in parts if p != rid]:
        if t and t not in seen:
            seen.add(t)
            out.append(t)
    return out


def canonical_gene_id(raw: str) -> str:
    """
    Strip a leading ``accession|`` prefix (text before the first ``|``).

    Wheat-style IDs without ``|`` are returned unchanged (aside from strip).
    """
    s = (raw or "").strip()
    if "|" not in s:
        return s
    _, _, rest = s.partition("|")
    rest = rest.strip()
    return rest if rest else s
