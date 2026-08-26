# PANDAGMA workflow

Cluster-ready wrappers around upstream [legumeinfo/pandagma](https://github.com/legumeinfo/pandagma) for building pan-gene clusters from annotated genomes.

## Installation

```bash
cd pandagma
make setup
conda activate pandagma
mkdir -p log
```

`make setup` clones `pandagma-upstream/`, applies local patches, and creates the conda environment from `environment.yml`.

Edit Slurm headers (`#SBATCH --account`, `#SBATCH --partition`) or pass overrides at submit time.

## Input layout

Create one directory per genome set (e.g. `wheat/`):

| Path | Content |
|------|---------|
| `wheat/cds/` | CDS FASTA per accession |
| `wheat/prot/` | Protein FASTA per accession |
| `wheat/config/wheat_pan.conf` | Pan config (+ optional `*_expected_chr_matches.inc`) |
| `wheat/work/` | Created by the pipeline |

BED files (seven columns) and pan config fields follow upstream Pandagma conventions.

## Build pan-gene clusters

```bash
sbatch slurm/pandagma_ingest.slurm wheat
sbatch slurm/pandagma_mmseqs_array.slurm wheat
bash scripts/submit_pandagma_dag_separate_jobs.sh wheat
sbatch slurm/pandagma_dagchainer_finalize.slurm wheat
sbatch slurm/pandagma_pan_resume.slurm wheat
```

Final output: `wheat/work/18_syn_pan_aug_extra.clust.tsv`

## Build a PanViewer database

From the PanViewer repo root:

```bash
bash scripts/stage_variant_db.sh wheat.pandagma \
  --pan-tsv pandagma/wheat/work/18_syn_pan_aug_extra.clust.tsv \
  --bed-dir pandagma/wheat/bed \
  --prot-dir pandagma/wheat/prot
```

This writes `database/wheat.pandagma.db`. Register the variant in `database/config.json` (see `database/config.json.example`).

## Scaling

Tune DAG submission with environment variables documented in `scripts/submit_pandagma_dag_separate_jobs.sh` (`DAG_MAX_JOBS`, `DAG_CPUS`, `DAG_MEM`, `DAG_TIME`).
