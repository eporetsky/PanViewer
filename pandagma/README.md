# PANDAGMA (Slurm workflow)

Cluster-ready wrappers around upstream [legumeinfo/pandagma](https://github.com/legumeinfo/pandagma) for building pan-gene clusters from a set of annotated genomes. The final tabular product used by PanViewer is a file like:

```text
<genome>/work/18_syn_pan_aug_extra.clust.tsv
```

Copy or symlink that TSV (plus `bed/`, optional `prot/` / `cds/`) into PanViewer’s `input/<species>/`, then run `python build_index.py --force` from the PanViewer repo root.

## Setup

From this directory (`pandagma/`):

```bash
make setup
conda activate pandagma
mkdir -p log
```

`make setup` clones upstream into `pandagma-upstream/`, applies local patches (ingest FASTA suffix + DAGchainer reliability), and creates/updates the `pandagma` conda env (`environment.yml`).

If `pandagma: command not found` after activate:

```bash
conda deactivate && conda activate pandagma
# or from this directory:
export PATH="$(pwd)/pandagma-upstream/bin:$PATH"
```

### Slurm account / partition

Scripts use placeholders. Edit them (or pass on the CLI) before submitting:

```bash
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=YOUR_PARTITION
```

Example:

```bash
sbatch --account=myacct --partition=myqueue slurm/pandagma_ingest.slurm wheat
```

## Directory layout (per genome set)

Create a directory named for your analysis (examples: `wheat`, `barley`, `my_crop`). Paths below assume `wheat/`:

| Path | Role |
|------|------|
| `wheat/cds/` | CDS FASTA per accession |
| `wheat/prot/` | Protein FASTA per accession |
| `wheat/config/` | Pan config: `wheat_pan.conf` or `pan.wheat.conf`, plus `*_expected_chr_matches.inc` when using chromosome filters |
| `wheat/work/` | Pandagma work dir (`02_fasta_nuc`, `03_mmseqs`, `04_dag`, …) |

Config discovery: `config/pan.<dirname>.conf` or `config/<dirname>_pan.conf`, or set `CONF=` to an absolute path. Work dir defaults to `<dirname>/work` (override with `WORK=`).

### Input conventions

- `annotation_files`, `cds_files`, and `protein_files` in the pan config must stay in **matching order**.
- BED: seven columns — chrom, start0, end, mRNA ID, score, strand, gene ID.
- Chromosome tokens used with `expected_chr_matches` should match BED molecules after stripping a leading `chr` (e.g. token `1A` ↔ `chr1A`).

## Wheat run (exact sequence)

This is the sequence that produced `wheat/work/18_syn_pan_aug_extra.clust.tsv`:

```bash
cd pandagma
mkdir -p log

sbatch slurm/pandagma_ingest.slurm wheat
# wait until ingest finishes

sbatch slurm/pandagma_mmseqs_array.slurm wheat
# wait until all mmseqs shard jobs finish
# (each shard also runs the chromosome-aware filter → 04_dag/*_matches.tsv)

bash slurm/submit_pandagma_dag_separate_jobs.sh wheat
# wait until all pdg-dag batch jobs finish

sbatch slurm/pandagma_dagchainer_finalize.slurm wheat
# builds 05_filtered_pairs.tsv

sbatch slurm/pandagma_pan_resume.slurm wheat
# mcl → tabularize (through add_extra); writes 18_syn_pan_aug_extra.clust.tsv under work/
```

Do **not** run `pandagma pan -s dagchainer` on a single node after the sharded DAG submit — that redoes all pairs.

## Custom genome set

1. **Prepare a directory** `pandagma/<name>/` with `cds/`, `prot/`, annotations/BED as required by your pan config, and `config/<name>_pan.conf` (or `pan.<name>.conf`).
2. **Tune the pan config** for your species (identity/coverage, `expected_chr_matches`, optional `dag_A` / `dag_M` / `dag_E` / `dag_gap_mult`, etc.). See upstream Pandagma docs for field meanings.
3. **Run the same five steps**, substituting your directory name for `wheat`:

```bash
sbatch slurm/pandagma_ingest.slurm <name>
sbatch slurm/pandagma_mmseqs_array.slurm <name>
bash slurm/submit_pandagma_dag_separate_jobs.sh <name>
sbatch slurm/pandagma_dagchainer_finalize.slurm <name>
sbatch slurm/pandagma_pan_resume.slurm <name>
```

4. **Feed PanViewer**: copy the resulting `*.clust.tsv` (and `bed/`, optional FASTA dirs) into `../input/<species>/`, then from the PanViewer root:

```bash
python build_index.py --force
```

### Scaling notes

- DAGchainer submit defaults (~48 CPU / 256GB per job, many separate jobs) suit medium–large pans. Tune with env vars on submit, e.g. `DAG_MAX_JOBS`, `DAG_PAIRS_PER_JOB`, `DAG_CPUS`, `DAG_MEM`, `DAG_TIME` (see header of `slurm/submit_pandagma_dag_separate_jobs.sh`).
- Optional: `PANDAGMA_CONDA_ENV=/path/to/conda/envs/pandagma` if cluster jobs need an explicit env prefix on `PATH`.

## What each Slurm step does

| Step | Script | Effect |
|------|--------|--------|
| Ingest | `slurm/pandagma_ingest.slurm` | `pandagma pan -s ingest` |
| Homology | `slurm/pandagma_mmseqs_array.slurm` | One mmseqs job per query genome; filter writes `04_dag/*_matches.tsv` |
| Synteny | `slurm/submit_pandagma_dag_separate_jobs.sh` | Many separate DAGchainer jobs over match files |
| Finalize | `slurm/pandagma_dagchainer_finalize.slurm` | Builds `05_filtered_pairs.tsv` |
| Resume | `slurm/pandagma_pan_resume.slurm` | `mcl` … `tabularize` (skips re-running dagchainer) |

## Local patches (applied by `make setup`)

- **Ingest FASTA suffix** — `scripts/patch-pandagma-pan-ingest-fasta.sh` so hashed CDS files end in `.fa` and downstream globs work.
- **DAGchainer batch reliability** — `scripts/patch-pandagma-pan-dagchainer.sh` for sharded DAGchainer runs.

Re-run `make patch-pandagma-ingest-fasta` / `make patch-pandagma-dagchainer` after refreshing `pandagma-upstream/`.
