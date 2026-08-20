# GeneTribe → pan-genes (PanViewer)

Pipeline from **protein FASTA + BED or GFF** to **pan-gene groups**, then an optional PanViewer SQLite annotation DB.

GeneTribe runs pairwise homology only. This repo adds the pan-gene step:

**Subgenome-aware + collinear:** keep an RBH only if both genes share a chromosome group across homeologous subgenomes (wheat `1A+1B+1D`, oat `1A+1C+1D`, barley `1H`…`7H`) **and** the pair sits in a collinear block. Then take connected components. That is the production rule used for `*.genetribe.db`.

---

## Setup

One conda env from `environment.yml` (blast, bedtools, jcvi, …). The `genetribe` CLI is **not** on conda — `make setup` clones it next to this folder and wires PATH.

```bash
cd genetribe
make setup                 # conda env + clone genetribe-upstream/
conda activate genetribe
mkdir -p log
which genetribe            # …/genetribe-upstream/genetribe
```

Optional overrides (export before running, or pass on the command line):

| Variable | Default / meaning |
|----------|-------------------|
| `PANVIEWER_ROOT` | Parent of `genetribe/` if it has `build_pangene_index.py` |
| `GT_CPUS` / `GT_MEM` / `GT_TIME` | Slurm pair-job resources (48 / 128G / 72h) |
| `SLURM_ACCOUNT` / `SLURM_PARTITION` | `small_grains` / `atlas` |
| `GENETRIBE_THREADS` | Local pair threads (default 8) |

---

## Directory layout (one species)

Example for `wheat/`:

| Path | What it is |
|------|------------|
| `wheat/prot/<acc>.fa` | Protein FASTA (gene-level headers) |
| `wheat/bed/<acc>.bed` **or** `wheat/gff/<acc>.gff3` | Coordinates |
| `wheat/accessions/<acc>/<acc>.{fa,bed,chrlist}` | Prepared GeneTribe inputs (from step 1) |
| `wheat/work/pairs/<A>_x_<B>/` | Pairwise GeneTribe outputs (`*.RBH`, collinearity, …) |
| `wheat/work/genetribe_subgenome_pans.hsh.tsv` | Pan membership (`pan_id` + `gene_id`) |
| `wheat/work/genetribe_subgenome_pans.clust.tsv` | Wide cluster table |

**Gene IDs:** FASTA headers and BED column 4 must be the **same gene-level ID**. Strip transcript suffixes yourself. This workflow always passes GeneTribe `-s @`.

**Chromosome names** in BED/GFF must look like `1A` / `chr1A` / `1H` (digit + subgenome letter).

---

## Step-by-step

### 1. Prepare accessions (protein + BED or GFF)

```bash
bash scripts/genetribe_prep.sh wheat \
  --chrlist examples/wheat.chrlist \
  --prot-dir wheat/prot --bed-dir wheat/bed
# or: --gff-dir wheat/gff

bash scripts/genetribe_check_ids.sh wheat
```

**Outputs:** `wheat/accessions/<acc>/<acc>.{fa,bed,chrlist}`

If chroms are `chr1A`… use a chrlist like `chrNA,chrNB,chrND` instead of `NA,NB,ND`.

### 2. Pairwise GeneTribe (`genetribe core`)

**Slurm:**

```bash
bash slurm/submit_pairs.sh wheat
# wait until: find wheat/work/pairs -name '*.RBH' | wc -l   ≈ N*(N-1)/2
```

**Local:**

```bash
GENETRIBE_THREADS=16 bash scripts/genetribe_run_pairs.sh wheat
```

**Per-pair outputs** under `work/pairs/<A>_x_<B>/`: `*.RBH` (used for pans), collinearity files, `.done`.

### 3. Cluster pan-genes (subgenome-aware + collinear)

```bash
sbatch slurm/cluster_pans.slurm wheat
# or: python3 scripts/genetribe_cluster_pans.py wheat
```

**Outputs:** `work/genetribe_subgenome_pans.{hsh,clust}.tsv`

### 4. Build PanViewer annotation DB (optional)

Needs `database/<species>.db` already built (`build_base_index.py`).

```bash
sbatch slurm/build_db.slurm wheat
# or: bash scripts/genetribe_build_db.sh wheat
# or: PANVIEWER_ROOT=/path/to/panviewer bash scripts/genetribe_build_db.sh wheat
```

**Output:** `database/<species>.genetribe.db`

---

## Script map

| Step | Local | Slurm |
|------|-------|-------|
| Setup | `make setup` | — |
| Prep | `scripts/genetribe_prep.sh` | — |
| ID check | `scripts/genetribe_check_ids.sh` | (also run by submit) |
| Pairs | `scripts/genetribe_run_pairs.sh` | `slurm/submit_pairs.sh` → `pair_job.slurm` |
| Cluster | `scripts/genetribe_cluster_pans.py` | `slurm/cluster_pans.slurm` |
| DB | `scripts/genetribe_build_db.sh` | `slurm/build_db.slurm` |

---

## Citations

- Chen et al. (2020) GeneTribe. *Molecular Plant*. https://chenym1.github.io/genetribe/
- Yu et al. (2023) Rice Gene Index. *Molecular Plant*. https://doi.org/10.1016/j.molp.2023.03.012
