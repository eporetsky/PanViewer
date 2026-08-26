# PanViewer

**PanViewer** is a web-based tool for exploring protein-level variation across a pangenome: find genes quickly, open any **pangene**, and move between interactive alignments, variant summaries, a neighbor-joining guide tree, and exports—without leaving the browser.

**Scope:** Pan-genes for **wheat**, **barley**, and **oat**. Each species can expose multiple **variants** (Pandagma, GeneTribe, OrthoFinder-derived sets) via `database/config.json`. Wheat additionally supports coordinate-aware synteny and cluster UI when `gene_coords` is present in the index.

## Features

- Gene search with autocomplete (maps each gene to its pangene), with cross-species ID hints and automatic dataset switching where patterns are recognizable
- **Keyword search** (Arabidopsis / rice annotation text) linked to pan-genes via mmseqs best hits, scoped to the active species tab
- **FAMSA** multiple sequence alignment with an interactive MSA viewer (reference switching, gap and color options)
- Variant summary vs. a chosen reference sequence
- Phylogenetic tree view (topology from the same FAMSA NJ guide tree)
- Downloads: gene lists, protein FASTA, alignment FASTA, variant TSV
- Saved genes / pangenes in the browser (`localStorage`)
- With **gene coordinates** in the index: local synteny context and subgenome-oriented gene tables (wheat A/B/D, oat A/C/D)
- **Expression tab** (wheat / barley / oat reference genes): tissue CPM + DEG plots from local `expression/expression.db` when present, otherwise [PlantApp](https://www.plantapp.org)

## Data layout

### Species databases (`database/`)

SQLite files built from staged inputs under `input/<staging_id>/`. Use **`database/config.json`** (copy from `database/config.json.example`) to group variants under species tabs:

| Staging id | DB file | Method (example) |
|------------|---------|------------------|
| `wheat.pandagma` | `wheat.pandagma.db` | Pandagma |
| `wheat.genetribe` | `wheat.genetribe.db` | GeneTribe |
| `barley.panbarlex` | `barley.panbarlex.db` | OrthoFinder (PanBARLEX) |
| `oat.panoat` | `oat.panoat.db` | OrthoFinder (PanOat) |

Each variant input tree holds a pan TSV (`*.clust.tsv` or `*.hsh.tsv`), a `bed/` directory, and optional `prot/`, `cds/`, `porter6/` (often symlinks to shared species-wide data). See **`database/README.md`** for customizing config.

| File | Role |
|------|------|
| `database/*.db` | One SQLite index per variant |
| `database/config.json` | Species tabs and variant labels (optional) |
| `database/stats.tsv` | One row per variant; updated by `build_index.py` or `update_database_stats.py` |

Without `config.json`, each `database/<stem>.db` becomes its own species tab.

### Keyword search (`search/`)

| Path | Role |
|------|------|
| `search/reference/` | Arabidopsis and rice annotation files (indexed once, shared across species) |
| `search/mmseqs_wheat/`, `search/mmseqs_barley/`, `search/mmseqs_oat/` | Per-accession TSVs: `gene_id`, `arabidopsis`, `rice`, optional e-value / pident columns |
| `search/keyword_index.sqlite` | **Single shared** FTS + `pan_hits` index (built by `build_keyword_index.py`) |

**One keyword index for all species** — not separate indexes per genome. `pan_hits` rows store `dataset_id` (`wheat`, `barley`, `oat`); `/search` filters keyword results to the active species tab. After adding or renaming mmseqs TSV directories, or rebuilding a species DB, **rebuild the keyword index** (see below).

Legacy layout `search/mmseqs/` (single directory) is still supported if present.

### Expression data (`expression/`)

Optional. If **`expression/expression.db`** exists, the Expression tab reads that file and **does not call** the PlantApp API. One DB can hold multiple datasets (wheat, barley, oat, or custom).

| Path | Role |
|------|------|
| `expression/<dataset>/meta.tsv` | Sample metadata (`sample_acc`, organ, group, titles, …) |
| `expression/<dataset>/cpm.tsv[.gz]` | Gene × sample CPM matrix |
| `expression/<dataset>/deg.tsv[.gz]` | Optional long DEG table |
| `expression/<dataset>/dataset.json` | Optional `genome` / `label` |
| `expression/expression.db` | Combined runtime DB (like `search/keyword_index.sqlite`; overrides PlantApp when present) |

See **`expression/README.md`** for column lists and size tips. Build with:

```bash
python build_expression.py --force
```

Without `expression.db`, PanViewer proxies PlantApp (`PLANTAPP_BASE_URL`, default `https://www.plantapp.org`).

## Installation

### Conda environment

```bash
conda create -n panviewer python=3.11 flask biopython famsa -c conda-forge -c bioconda -y
conda activate panviewer
```

### Build PANDAGMA pan-genes (optional; custom genomes)

To create the pan-gene cluster TSV from scratch on a Slurm cluster (wheat and other custom genome sets), see **[`pandagma/README.md`](pandagma/README.md)**. Short form:

```bash
cd pandagma
make setup && conda activate pandagma && mkdir -p log
# prepare <name>/{cds,prot,config}/ then:
sbatch slurm/pandagma_ingest.slurm <name>
sbatch slurm/pandagma_mmseqs_array.slurm <name>
bash scripts/submit_pandagma_dag_separate_jobs.sh <name>
sbatch slurm/pandagma_dagchainer_finalize.slurm <name>
sbatch slurm/pandagma_pan_resume.slurm <name>
```

Copy the resulting `<name>/work/18_syn_pan_aug_extra.clust.tsv` (and BED/FASTA inputs) into `input/<species>/` as below.

### Build GeneTribe pan-genes (optional; RGI-style)

Alternative to Pandagma: all-vs-all [GeneTribe](https://chenym1.github.io/genetribe/) pairs, then RBH connected-component clustering (same approach as Rice Gene Index). See **[`genetribe/README.md`](genetribe/README.md)**. Short form:

```bash
cd genetribe
make setup && conda activate genetribe && mkdir -p log
# prepare <name>/accessions/<acc>/<acc>.{fa,bed,chrlist} then:
bash slurm/submit_genetribe_pairs.sh <name>
sbatch slurm/finalize.slurm <name>
```

Copy `<name>/work/genetribe_pans.hsh.tsv` (or `.clust.tsv`) plus BED/FASTA into a variant staging tree, then build the database (step below).

### OrthoFinder pan-genes (optional)

If you run [OrthoFinder](https://github.com/davidemms/OrthoFinder) separately, convert HOG tables with **`orthofinder/reformat_orthofinder_to_pandagma.py`** — see **[`orthofinder/README.md`](orthofinder/README.md)**.

### Porter6 localization (optional)

Preprocess and run [Porter6](https://github.com/WafaAlanazi/Porter6) predictions for inclusion in variant indexes — see **[`porter6/README.md`](porter6/README.md)**.

### Build variant indexes

Stage method outputs and build SQLite indexes with the shared helper:

```bash
cp database/config.json.example database/config.json   # optional; customize variant labels

bash scripts/stage_variant_db.sh wheat.pandagma \
  --pan-tsv pandagma/wheat/work/18_syn_pan_aug_extra.clust.tsv \
  --bed-dir /path/to/wheat/bed \
  --prot-dir /path/to/wheat/prot

bash scripts/stage_variant_db.sh wheat.genetribe \
  --pan-tsv genetribe/wheat/work/genetribe_pans.hsh.tsv \
  --bed-dir genetribe/wheat/bed \
  --prot-dir genetribe/wheat/prot
```

Or call method-specific wrappers (e.g. `genetribe/scripts/genetribe_build_variant_db.sh wheat`). This writes `database/<staging_id>.db` and refreshes `database/stats.tsv`.

To rebuild one variant only:

```bash
python build_index.py --force wheat.pandagma
```

After a PanViewer upgrade that changes the SQLite layout, **rebuild every variant database** with `--force`.

### Build keyword search index

Requires variant DBs for species you want keyword search on, plus mmseqs TSVs under `search/mmseqs_*`:

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

### Build expression database (optional)

Place curated TSV/CSV folders under `expression/<dataset>/`, then:

```bash
python build_expression.py --force
```

This writes `expression/expression.db`. When that file is present, expression requests use it instead of PlantApp. Details: `expression/README.md`.

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

Place `database/*.db` (and optionally `search/keyword_index.sqlite`, `expression/expression.db`) before the build, or mount them at runtime. Behind a reverse proxy, set **`SCRIPT_NAME`** or **`X-Forwarded-Prefix`** to match the public path so generated links stay correct.

## Usage

1. Choose species (wheat / barley / oat) if more than one DB is present.
2. Search for a gene ID (partial IDs allowed where supported) or a keyword (gene name, GO text, etc.).
3. Open a pangene from the results.
4. Use the alignment, variants, tree, and download tabs as needed; bookmark items from the UI for quick return visits.

## Project layout

```
PanViewer/
├── app.py                 # Flask app
├── build_index.py         # SQLite DBs from input/<staging_id>/ → database/<staging_id>.db
├── build_keyword_index.py # Shared search/keyword_index.sqlite from reference + mmseqs_*
├── build_expression.py    # expression/<dataset>/ → expression/expression.db
├── scripts/stage_variant_db.sh  # Stage method outputs + build one variant DB
├── expression_local.py    # Local omics lookup (overrides PlantApp when DB exists)
├── keyword_search.py      # Keyword fallback for /search
├── plantapp_omics.py      # Local DB or PlantApp gene-omics proxy
├── dataset_stats.py
├── database_config.py     # Load database/config.json
├── environment.yml
├── pandagma/              # Slurm PANDAGMA workflow (see pandagma/README.md)
├── genetribe/             # Slurm GeneTribe → pan-genes (see genetribe/README.md)
├── orthofinder/           # OrthoFinder HOG → pan TSV (see orthofinder/README.md)
├── porter6/               # Porter6 batch helpers (see porter6/README.md)
├── database/              # *.db, config.json, stats.tsv
├── expression/            # per-dataset sources + expression.db
├── input/<staging_id>/    # Staged pan TSV + bed/prot/cds/porter6 per variant
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
