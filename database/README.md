# Database configuration

PanViewer discovers SQLite files under `database/` and optionally groups them into species tabs via `config.json`.

## Quick start

1. Copy the example config:

   ```bash
   cp database/config.json.example database/config.json
   ```

2. Build variant databases with `scripts/stage_variant_db.sh` or method-specific helpers (see `pandagma/`, `genetribe/`, `orthofinder/` READMEs).

3. Ensure each variant’s `db` filename matches what you built (e.g. `wheat.pandagma.db`).

## Customize config

Copy `database/config.json.example` to `database/config.json`. Each species tab lists **variants** (pangene methods). Two ids are involved:

| Role | Example | Meaning |
|------|---------|---------|
| **Staging id** | `wheat.pandagma` | `input/<id>/` folder and `database/<id>.db` filename (used by `scripts/stage_variant_db.sh`) |
| **Config variant key** | `wheat` | UI selector id; `"db"` must match the staged filename (`wheat.pandagma.db`) |

Edit labels, add or remove variants, and set `keyword_dataset_id` when keyword search should use a different species’ mmseqs index.

If `config.json` is absent, each `database/<stem>.db` becomes its own species tab.

## Stats

`database/stats.tsv` is updated by `build_index.py`. With `config.json` present, stats rows use variant ids from the config.
