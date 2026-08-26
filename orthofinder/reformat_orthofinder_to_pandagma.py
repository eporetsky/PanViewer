#!/usr/bin/env python3
"""Convert OrthoFinder HOG TSV to pandagma-style cluster TSV for ``build_index.py``.

Accepts:
  - raw OrthoFinder ``N0.tsv`` (columns: HOG, OG, Gene Tree Parent Clade, species…)
  - already-normalized species **table** (HOG + one column per accession)

``clust`` (default) is what PanViewer needs: pan-id, then one gene ID per tab field.
``table`` keeps species columns (the format previously uploaded as
``18_syn_pan_aug_extra.clust.tsv`` — that is *not* a pandagma clust file).

Gene cells may contain comma-separated paralogs. Transcript suffixes (``.1``) are
stripped. Species names ending in ``.primary.protein`` are shortened.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

TRANSCRIPT_SUFFIX_RE = re.compile(r"\.\d+$")
GENOME_NAME_SUFFIX = ".primary.protein"
RAW_N0_META = {"OG", "Orthogroup"}


def normalize_genome_name(name: str) -> str:
    if name.endswith(GENOME_NAME_SUFFIX):
        return name[: -len(GENOME_NAME_SUFFIX)]
    return name


def strip_transcript_id(gene_id: str) -> str:
    return TRANSCRIPT_SUFFIX_RE.sub("", gene_id.strip())


def parse_gene_cell(cell: str) -> list[str]:
    if not cell.strip():
        return []
    return [gene.strip() for gene in cell.split(",") if gene.strip()]


def species_column_start(header: list[str]) -> int:
    """Return first species-column index: 3 for raw N0, 1 for HOG+species table."""
    if not header or header[0] != "HOG":
        raise ValueError(
            "Unexpected input format: first column must be HOG "
            "(OrthoFinder N0.tsv or a HOG+species table)."
        )
    if len(header) >= 4 and header[1] in RAW_N0_META:
        return 3
    if len(header) < 2:
        raise ValueError("HOG table has no species columns.")
    return 1


def collect_hog_genes(row: list[str], species_start: int) -> list[str]:
    seen: set[str] = set()
    genes: list[str] = []
    for cell in row[species_start:]:
        for raw_gene in parse_gene_cell(cell):
            gene = strip_transcript_id(raw_gene)
            if gene not in seen:
                seen.add(gene)
                genes.append(gene)
    return genes


def format_species_cell(cell: str) -> str:
    seen: set[str] = set()
    genes: list[str] = []
    for raw_gene in parse_gene_cell(cell):
        gene = strip_transcript_id(raw_gene)
        if gene not in seen:
            seen.add(gene)
            genes.append(gene)
    return ",".join(genes)


def convert_orthofinder_to_pandagma(
    input_path: Path,
    output_path: Path,
    *,
    output_format: str = "clust",
    keep_singletons: bool = False,
) -> dict[str, int]:
    stats = {
        "input_hogs": 0,
        "skipped_single_gene": 0,
        "output_hogs": 0,
        "species_start": 0,
    }

    with input_path.open(newline="") as infile, output_path.open("w", newline="") as outfile:
        reader = csv.reader(infile, delimiter="\t")
        writer = csv.writer(outfile, delimiter="\t", lineterminator="\n")

        header = next(reader)
        species_start = species_column_start(header)
        stats["species_start"] = species_start
        species_names = [normalize_genome_name(col) for col in header[species_start:]]
        if output_format == "table":
            writer.writerow(["HOG", *species_names])

        min_genes = 1 if keep_singletons else 2
        for row in reader:
            if not row or not (row[0] or "").strip():
                continue
            stats["input_hogs"] += 1
            hog_id = row[0].strip()
            genes = collect_hog_genes(row, species_start)
            if len(genes) < min_genes:
                stats["skipped_single_gene"] += 1
                continue
            if output_format == "table":
                cells = row[species_start:]
                if len(cells) < len(species_names):
                    cells = cells + [""] * (len(species_names) - len(cells))
                writer.writerow([hog_id, *[format_species_cell(c) for c in cells]])
            else:
                writer.writerow([hog_id, *genes])
            stats["output_hogs"] += 1

    return stats


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Convert OrthoFinder HOG TSV (raw N0 or species table) to pandagma clust TSV "
            "for build_index.py."
        )
    )
    parser.add_argument(
        "input",
        type=Path,
        help="OrthoFinder N0.tsv or HOG+species table TSV",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="Output TSV (use *.clust.tsv for build_index.py)",
    )
    parser.add_argument(
        "-f",
        "--format",
        choices=("clust", "table"),
        default="clust",
        help="clust: pan-id + one gene per field (default). table: HOG + species columns.",
    )
    parser.add_argument(
        "--keep-singletons",
        action="store_true",
        help="Keep HOGs with a single gene (default: skip them).",
    )
    args = parser.parse_args()

    if not args.input.is_file():
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        return 1

    stats = convert_orthofinder_to_pandagma(
        args.input,
        args.output,
        output_format=args.format,
        keep_singletons=args.keep_singletons,
    )
    kind = "raw N0" if stats["species_start"] == 3 else "HOG+species table"
    print(f"Detected {kind} (species columns start at {stats['species_start']})")
    print(f"Wrote {stats['output_hogs']:,} clusters to {args.output}")
    print(
        f"Skipped {stats['skipped_single_gene']:,} HOGs with "
        f"<{'1' if args.keep_singletons else '2'} genes "
        f"out of {stats['input_hogs']:,} input HOGs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
