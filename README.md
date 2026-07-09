# PanViewer

**PanViewer** is a web-based tool for exploring protein-level variation across a pangenome: find genes quickly, open any **pangene**, and move between interactive alignments, variant summaries, a neighbor-joining guide tree, and exports—without leaving the browser.

**Scope:** PANDAGMA pan-genes for **wheat**, **barley**, and **oat** (one SQLite database per species). Wheat additionally supports coordinate-aware synteny and cluster UI when `gene_coords` is present in the index.

## Features

- Gene search with autocomplete (maps each gene to its pangene), with cross-species ID hints and automatic dataset switching where patterns are recognizable
- **Keyword search** (Arabidopsis / rice annotation text) linked to pan-genes via mmseqs best hits, scoped to the active species tab
- **FAMSA** multiple sequence alignment with an interactive MSA viewer (reference switching, gap and color options)
- Variant summary vs. a chosen reference sequence
- Phylogenetic tree view (topology from the same FAMSA NJ guide tree)
- Downloads: gene lists, protein FASTA, alignment FASTA, variant TSV
- Saved genes / pangenes in the browser (`localStorage`)
- With **gene coordinates** in the index: local synteny context and subgenome-oriented gene tables (wheat A/B/D, oat A/C/D)

## Data layout

### Species databases (`database/`)

Runtime SQLite files (built from `input/<species>/`):

| File | Source |
|------|--------|
| `database/wheat.db` | `input/wheat/` (PANDAGMA clust TSV + `bed/`, optional `prot/`, `cds/`, `porter6/`) |
| `database/barley.db` | `input/barley/` |
| `database/oat.db` | `input/oat/` |
| `database/stats.tsv` | One row per species (accessions, pan-genes, genes, optional table counts); updated by `build_index.py` or `update_database_stats.py` |

The app discovers any `database/<species>.db` at startup (no server restart needed after adding a new file).

### Keyword search (`search/`)

| Path | Role |
|------|------|
| `search/reference/` | Arabidopsis and rice annotation files (indexed once, shared across species) |
| `search/mmseqs_wheat/`, `search/mmseqs_barley/`, `search/mmseqs_oat/` | Per-accession TSVs: `gene_id`, `arabidopsis`, `rice`, optional e-value / pident columns |
| `search/keyword_index.sqlite` | **Single shared** FTS + `pan_hits` index (built by `build_keyword_index.py`) |

**One keyword index for all species** — not separate indexes per genome. `pan_hits` rows store `dataset_id` (`wheat`, `barley`, `oat`); `/search` filters keyword results to the active species tab. After adding or renaming mmseqs TSV directories, or rebuilding a species DB, **rebuild the keyword index** (see below).

Legacy layout `search/mmseqs/` (single directory) is still supported if present.

## Installation

### Conda environment

```bash
conda create -n panviewer python=3.11 flask biopython famsa -c conda-forge -c bioconda -y
conda activate panviewer
```

### Build species indexes (PANDAGMA)

For each species under `input/<species>/`, place a Pandagma pan TSV (`*.clust.tsv` or `*.hsh.tsv`), a `bed/` directory, and optional `prot/`, `cds/`, `porter6/` FASTA dirs. Then:

```bash
python build_index.py --force
```

This writes `database/<species>.db` for every valid `input/<species>/` tree and refreshes `database/stats.tsv`.

To build one species only, keep a single subdirectory under `input/` or point inputs via the helpers in `build_index.py` (see script and `AGENTS.md`).

The SQLite layout stores protein, CDS, and Porter6 in deduplicated **`*_uniq` / `*_map`** table pairs. After a PanViewer upgrade that changes this layout, **rebuild every species database** with `--force`.

For fast synteny-style window queries on wheat (or any species with coordinates), the database needs **`gene_coords.chrom_index`** and the window index on `(accession, chr, chrom_index)` — included in a normal PANDAGMA build with BED inputs.

### Build keyword search index

Requires `database/*.db` for species you want keyword search on, plus mmseqs TSVs under `search/mmseqs_*`:

```bash
python build_keyword_index.py --force
```

By default this ingests **all** `search/mmseqs_*` directories and matches each `<accession>.tsv` to the species DB that contains that accession. Override with repeated `--mmseqs-dir` flags, e.g.:

```bash
python build_keyword_index.py --force \
  --mmseqs-dir search/mmseqs_wheat \
  --mmseqs-dir search/mmseqs_barley \
  --mmseqs-dir search/mmseqs_oat
```

To regenerate mmseqs TSVs from protein FASTA (optional; slow), see `search/map_prot_to_at_os.py`.

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

Place `database/*.db` (and optionally `search/keyword_index.sqlite`) before the build, or mount them at runtime. Behind a reverse proxy, set **`SCRIPT_NAME`** or **`X-Forwarded-Prefix`** to match the public path so generated links stay correct.

## Usage

1. Choose species (wheat / barley / oat) if more than one DB is present.
2. Search for a gene ID (partial IDs allowed where supported) or a keyword (gene name, GO text, etc.).
3. Open a pangene from the results.
4. Use the alignment, variants, tree, and download tabs as needed; bookmark items from the UI for quick return visits.

## Project layout

```
PanViewer/
├── app.py                 # Flask app
├── build_index.py         # SQLite DBs from input/<species>/ → database/<species>.db
├── build_keyword_index.py # Shared search/keyword_index.sqlite from reference + mmseqs_*
├── keyword_search.py      # Keyword fallback for /search
├── dataset_stats.py
├── environment.yml
├── database/              # wheat.db, barley.db, oat.db, stats.tsv
├── input/<species>/       # PANDAGMA sources per species
├── search/
│   ├── reference/
│   ├── mmseqs_wheat/ | mmseqs_barley/ | mmseqs_oat/
│   └── keyword_index.sqlite
├── templates/
└── static/
```

## Notes

- **Alignment:** `famsa` on `PATH` (FAMSA2: neighbour-joining guide via ``-gt nj -gt_export`` then ``-gt import``; displayed tree is that NJ Newick **midpoint-rooted** in Python). If NJ export/import is unavailable, the app falls back to default FAMSA plus **FastTree** on the MSA (also midpoint-rooted).
- **Consensus:** Majority rule over non-gap residues at each column.
- **Presets / bookmarks:** Stored only in the browser; nothing is written on the server for those features.
