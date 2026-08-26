# Porter6 preprocessing

Helpers for running [Porter6](https://github.com/WafaAlanazi/Porter6) subcellular-localization predictions on large pangenome protein sets and expanding results back to per-genotype CSV files for PanViewer.

Porter6 itself is **not** vendored here — install it separately and point `PORTER6_ROOT` at your checkout.

## Installation

1. Clone Porter6 and install its dependencies (CUDA, PyTorch, etc.).
2. Place or symlink the Porter6 code under `porter6/porter6/` in this directory.
3. Create a conda env with Biopython and pandas for the preprocessing scripts.

## Workflow

### 1. Deduplicate sequences (optional, for pangenomes)

```bash
python scripts/pangenome_dedupe.py \
  --input-dir /path/to/prot \
  --out-dir pangenome/
```

### 2. Split and run Porter6 on a cluster

```bash
python scripts/split_fasta_chunks.py \
  --input pangenome/unique.fasta \
  --out-dir batches/ --chunk-size 100000

bash slurm/porter6.batch.sh batches/
```

Each job writes `results/<stem>/porter6_predictions_merged.csv`. Intermediate `test_ensemble_*.json` files can be deleted after recombination.

### 3. Expand to per-genotype CSVs

```bash
python scripts/recombine_predictions.py \
  --map-tsv pangenome/orig_to_uniq.tsv \
  --predictions-dir results/ \
  --out-dir porter6/
```

Outputs `{genotype}.q3.csv` and `{genotype}.q8.csv`.

### 4. Add to PanViewer

Symlink or copy the CSV directory into your variant input tree:

```bash
# when staging a DB variant:
bash scripts/stage_variant_db.sh wheat.pandagma \
  --pan-tsv ... --bed-dir ... --prot-dir ... \
  --porter6-dir porter6/
```

## Slurm

Edit `#SBATCH --account` and `#SBATCH --partition` in `slurm/porter6.sh`. Set `PORTER6_CUDA_MODULE` (e.g. `cuda/12.4.0`) if your site uses environment modules. GPU partition and memory requirements depend on your Porter6 install.
