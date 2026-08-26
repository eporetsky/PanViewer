# OrthoFinder → PanViewer

Convert OrthoFinder HOG tables into a pandagma-style cluster TSV for `build_index.py`.

This repo does **not** include OrthoFinder itself — run OrthoFinder separately, then use the script below.

## Usage

```bash
python orthofinder/reformat_orthofinder_to_pandagma.py \
  Orthogroups/N0.tsv \
  -o input/oat.panoat/panoat.clust.tsv \
  --format clust
```

Accepts raw OrthoFinder `N0.tsv` or an already-normalized HOG+species table. Singleton HOGs are skipped by default; use `--keep-singletons` to retain them.

## Build PanViewer database

```bash
bash scripts/stage_variant_db.sh oat.panoat \
  --pan-tsv input/oat.panoat/panoat.clust.tsv \
  --bed-dir /path/to/oat/bed \
  --prot-dir /path/to/oat/prot
```

Writes `database/oat.panoat.db`. Register variants in `database/config.json` (see `database/config.json.example`).

For barley OrthoFinder HOGs (PanBARLEX), use staging id `barley.panbarlex` → `database/barley.panbarlex.db`.
