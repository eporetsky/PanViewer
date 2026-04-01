# PanViewer

**PanViewer** is a web-based tool for exploring protein-level variation across a pangenome: find genes quickly, open any **pangene**, and move between interactive alignments, variant summaries, a neighbor-joining guide tree, and exports—without leaving the browser.

**Scope:** The pipeline and UI are currently oriented toward **hexaploid wheat** (pangenes prepared with [PANDAGMA](https://github.com/legumeinfo/pandagma)). Later releases aim to make the same workflow usable for other species.

## Features

- Gene search with autocomplete (maps each gene to its pangene)
- **FAMSA** multiple sequence alignment with an interactive MSA viewer (reference switching, gap and color options)
- Variant summary vs. a chosen reference sequence
- Phylogenetic tree view (topology from the same FAMSA NJ guide tree)
- Downloads: gene lists, protein FASTA, alignment FASTA, variant TSV
- Saved genes / pangenes in the browser (`localStorage`)
- With **gene coordinates** in the index: local synteny context and wheat subgenome (A/B/D)–oriented tables where applicable

## Data layout (`input/`)

- **`input/wheat/panwheat_pandagma.db`** — SQLite index produced by `build_index.py` from PANDAGMA pan membership and BED coordinates (see below). The app expects this path unless you change `DATASETS` in `app.py`.
- **`input/wheat/dataset_stats.tsv`** — optional one-row summary for the home page (accessions, pangenes, genes). Written on a full `build_index.py` run when applicable.

For fast synteny-style window queries, the database should include **`gene_coords.chrom_index`** and the window index on `(accession, chr, chrom_index)`. A normal PANDAGMA build with `--pandagma-bed-dir` handles this; to refresh indexes on an existing DB:

```bash
python build_index.py --chrom-index-only --out input/wheat/panwheat_pandagma.db
```

## Installation

### Conda environment

```bash
conda create -n panviewer python=3.11 flask biopython famsa -c conda-forge -c bioconda -y
conda activate panviewer
```

### Build the index (PANDAGMA)

Point at your pan membership TSV (wide clust or two-column `pan_id` / `gene_id`) and a directory of per-accession BED files. Protein (and optionally CDS) FASTA directories default under `primary/` if present; override with `--pandagma-prot-dir` / `--pandagma-cds-dir`.

```bash
python build_index.py \
  --pandagma-pan-tsv path/to/clust.tsv \
  --pandagma-bed-dir path/to/bed \
  --out input/wheat/panwheat_pandagma.db \
  --force
```

If your PANDAGMA outputs already live at the paths expected by this repo (e.g. `primary/…`), you can also run `python build_index.py` or `python build_index.py --force` with no extra flags and inspect the script’s console output.

### Run locally

```bash
python app.py
```

Open the app using the URL prefix set in `app.py` as **`APPLICATION_ROOT`** (e.g. `http://localhost:5050/panviewer/` if unchanged).

## Docker

From the repository root:

```bash
docker build -t panviewer .
docker run -d --name panviewer -p 8080:80 --restart unless-stopped panviewer
```

Place `input/wheat/` (including the SQLite DB) before the build, or mount it at runtime. Behind a reverse proxy, set **`SCRIPT_NAME`** or **`X-Forwarded-Prefix`** to match the public path so generated links stay correct.

## Usage

1. Search for a gene ID (partial IDs allowed where supported).
2. Open a pangene from the results.
3. Use the alignment, variants, tree, and download tabs as needed; bookmark items from the UI for quick return visits.

## Project layout

```
PanViewer/
├── app.py                 # Flask app
├── build_index.py         # SQLite index from PANDAGMA inputs (+ optional upgrades)
├── dataset_stats.py       # Stats helpers
├── environment.yml
├── input/wheat/           # Wheat DB and optional dataset_stats.tsv
├── templates/
└── static/
```

## Notes

- **Alignment:** `famsa` must be on `PATH`. The server uses FAMSA’s NJ guide-tree export/import so the MSA and tree stay consistent.
- **Consensus:** Majority rule over non-gap residues at each column.
- **Presets / bookmarks:** Stored only in the browser; nothing is written on the server for those features.
