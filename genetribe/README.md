# GeneTribe workflow

Pipeline from protein FASTA + BED/GFF to pan-gene groups, then a PanViewer SQLite variant database.

GeneTribe performs pairwise homology; this repo adds RBH connected-component clustering (`genetribe_rbh_to_pans.py`).

## Installation

```bash
cd genetribe
make setup
conda activate genetribe
mkdir -p log
```

## Input layout

Example for `wheat/`:

| Path | Content |
|------|---------|
| `wheat/prot/<acc>.fa` | Protein FASTA (gene-level headers) |
| `wheat/bed/<acc>.bed` or `wheat/gff/` | Coordinates |
| `wheat/accessions/<acc>/` | Prepared GeneTribe inputs (step 1) |
| `wheat/work/pairs/` | Pairwise GeneTribe outputs |

Gene IDs must match between FASTA and BED. Use six-column BED (chrom, start, end, gene_id, score, strand).

## Step 1 — Prepare accessions

```bash
bash scripts/genetribe_prep.sh wheat \
  --chrlist examples/wheat.chrlist \
  --prot-dir wheat/prot --bed-dir wheat/bed
bash scripts/genetribe_check_ids.sh wheat
```

## Step 2 — Pairwise GeneTribe

```bash
bash slurm/submit_pairs.sh wheat
# or locally: bash scripts/genetribe_run_pairs.sh wheat
```

## Step 3 — Cluster pan-genes

```bash
sbatch slurm/finalize.slurm wheat
```

Writes `wheat/work/genetribe_pans.{hsh,clust}.tsv`.

## Step 4 — Build PanViewer database

From the PanViewer repo root:

```bash
bash genetribe/scripts/genetribe_build_variant_db.sh wheat
```

Writes `database/wheat.genetribe.db` (staging id `wheat.genetribe`). Add the variant to `database/config.json` (see `database/config.json.example`).

## Alternative clustering

`scripts/genetribe_cluster_pans.py` applies subgenome-aware + collinear filtering (see `slurm/cluster_pans.slurm`). Outputs `genetribe_subgenome_pans.*` instead of `genetribe_pans.*`.

## Citations

- Chen et al. (2020) GeneTribe. https://chenym1.github.io/genetribe/
- Yu et al. (2023) Rice Gene Index. https://doi.org/10.1016/j.molp.2023.03.012
