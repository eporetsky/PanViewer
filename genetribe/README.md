# GeneTribe → PanViewer pan-genes

Slurm-oriented workflow that runs [GeneTribe](https://chenym1.github.io/genetribe/) pairwise homology, then builds pan-gene groups with the Rice Gene Index strategy ([Yu et al. 2023](https://doi.org/10.1016/j.molp.2023.03.012)): **connected components on RBH edges only**.

1. Precompute directed BLASTPs into `work/blast/` (recommended).
2. Run `genetribe core -d work/blast` for every unordered accession pair.
3. Cluster RBH edges into pan-genes; write PanViewer-ready TSVs.

GeneTribe does not emit pan-gene clusters by itself — the finalize step is required.

## Requirements

| Component | Role |
|-----------|------|
| Conda env from `environment.yml` | `blast`, `bedtools`, `jcvi` (MCscan), Python |
| `make setup` → `genetribe-upstream/` | GeneTribe CLI (`install.sh`) — **not** a conda package |

```bash
make setup
conda activate genetribe
mkdir -p log
which genetribe   # …/genetribe-upstream/genetribe
```

## Analysis layout

Per analysis directory (example `wheat/`):

| Path | Role |
|------|------|
| `accessions/<acc>/<acc>.fa` | Protein FASTA |
| `accessions/<acc>/<acc>.bed` | 6-column gene BED |
| `accessions/<acc>/<acc>.chrlist` | Chromosome-group patterns ([formats](https://chenym1.github.io/genetribe/tutorial/fileformats.html)) |
| `work/blast/` | Precomputed `A_B.blast` files + BLAST DBs |
| `work/pairs/<A>_x_<B>/` | Per-pair GeneTribe outputs |
| `work/genetribe_pans.hsh.tsv` | PanViewer membership (`pan_id`, `gene_id`) |
| `work/genetribe_pans.clust.tsv` | Wide cluster table |

### Prep from Pandagma-style `prot/` + `bed/`

```bash
# chrlist uses GeneTribe's N placeholder, e.g. wheat:
#   echo 'chrNA,chrNB,chrND' > wheat/wheat.chrlist
bash scripts/genetribe_prep_from_pandagma.sh wheat --chrlist wheat/wheat.chrlist
```

Prep builds `accessions/` and converts 7-column Pandagma BEDs to GeneTribe’s 6 columns (`chrom start end gene_id score strand`). If `accessions/` already exists with 7-column beds:

```bash
bash scripts/genetribe_fix_beds.sh wheat
```

### Gene IDs

FASTA headers and BED column 4 must be the **same gene-level ID**. This workflow defaults to **`-s @`** so GeneTribe does not strip `.` inside real gene names (its default `-s .` will). Override only if needed:

```bash
GENETRIBE_EXTRA='-s |' bash slurm/submit_genetribe_pairs.sh wheat
```

Preflight (also run by submitters):

```bash
bash scripts/genetribe_check_id_separator.sh wheat
```

## Why precompute BLAST

Without a shared store, each pair job runs four BLASTPs (`A→B`, `B→A`, `A→A`, `B→B`). Self-blasts are repeated for every pair involving that accession. GeneTribe’s `core -d <dir>` reuses:

```text
work/blast/A_B.blast
work/blast/B_A.blast
work/blast/A_A.blast
work/blast/B_B.blast
```

This repo precomputes the full directed matrix (**N²** files, including self). Pair jobs then skip BLAST when all four files exist.

| Stage | Job count | What each job does |
|-------|----------:|--------------------|
| BLAST (recommended) | **N²** | One directed `blastp` |
| Pairs | **N(N−1)/2** | `genetribe core` (scoring + collinearity) |

`longestfasta` runs when building the blast store with the same `-s` as pair jobs (default `@`). With gene-level IDs this is usually a no-op; it keeps BLAST IDs consistent with GeneTribe.

## Slurm

Site-specific account/partition are **not** hardcoded. Set them if your cluster requires them:

```bash
export GT_ACCOUNT=myaccount
export GT_PARTITION=mypartition
# optional:
export GENETRIBE_CONDA_ENV=/path/to/conda/envs/genetribe
```

| Step | Script | Default CPUs / mem / time | Notes |
|------|--------|---------------------------|-------|
| BLAST | `slurm/submit_genetribe_blasts.sh` | 24 / 64G / 12h | One directed BLASTP per job |
| Pairs (blast store complete) | `slurm/submit_genetribe_pairs.sh` | 16 / 64G / 12h | Auto when `check_blasts` passes |
| Pairs (no/incomplete store) | same | 48 / 128G / 72h | Recomputes missing BLASTPs |
| Finalize | `slurm/genetribe_finalize.slurm` | 4 / 32G / 4h | After all `*.RBH` exist |

Override with `GT_CPUS`, `GT_MEM`, `GT_TIME`.

### Recommended run

```bash
mkdir -p log

bash slurm/submit_genetribe_blasts.sh wheat
# wait until the queue drains, then:
bash scripts/genetribe_check_blasts.sh wheat

bash slurm/submit_genetribe_pairs.sh wheat
# wait until:
#   find wheat/work/pairs -name '*.RBH' | wc -l   # ≈ N*(N-1)/2

sbatch slurm/genetribe_finalize.slurm wheat
# or: sbatch --account=... --partition=... slurm/genetribe_finalize.slurm wheat
```

Both submitters default to `GENETRIBE_SKIP_EXISTING=1` (skip finished BLAST files / completed pairs). Force reruns with `GENETRIBE_SKIP_EXISTING=0`.

### Monitoring

```bash
squeue -u $USER -h | wc -l
find wheat/work/blast -name '*.blast' -type f -size +0 | wc -l
bash scripts/genetribe_check_blasts.sh wheat

sstat -j <JOBID>.batch --format=JobID,AveCPU,MaxRSS,AveRSS,NTasks
srun --jobid=<JOBID> --overlap top -b -n1 | head -30
```

BLASTP often uses most but not all allocated cores (e.g. ~1800% on a 24-CPU job); that is normal.

## Hand off to PanViewer

```bash
# from the PanViewer repository
mkdir -p input/wheat_gt/bed
cp path/to/wheat/work/genetribe_pans.hsh.tsv input/wheat_gt/
# also copy bed files into input/wheat_gt/bed/
python build_index.py --force
```

## Pipeline scripts

| Step | Script |
|------|--------|
| Prep | `scripts/genetribe_prep_from_pandagma.sh` |
| Fix 7→6 col beds | `scripts/genetribe_fix_beds.sh` |
| BLAST store | `slurm/submit_genetribe_blasts.sh` |
| Check BLAST | `scripts/genetribe_check_blasts.sh` |
| Pair homology | `slurm/submit_genetribe_pairs.sh` |
| Cluster pans | `slurm/genetribe_finalize.slurm` |

Genes never seen in any RBH edge become size-1 pan-genes (from BED IDs). SBH / one2many tables under `work/pairs/` are kept for inspection but are **not** used for pan membership.

## Citations

- Chen et al. (2020) GeneTribe. *Molecular Plant*. https://chenym1.github.io/genetribe/
- Yu et al. (2023) Rice Gene Index. *Molecular Plant*. https://doi.org/10.1016/j.molp.2023.03.012
