#!/usr/bin/env python3
"""Expand unique-sequence predictions back to per-genome Q3/Q8 CSV files.

Outputs:
  recombined/<genome_id>.q3.csv   columns: id,q3_decoded
  recombined/<genome_id>.q8.csv   columns: id,q8_decoded
"""
from __future__ import annotations

import argparse
from glob import glob
from pathlib import Path

import pandas as pd


def genome_id_from_source(source_file: str) -> str:
    name = Path(str(source_file)).name
    for ext in (".fasta", ".faa", ".fa", ".FASTA", ".FAA", ".FA"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def load_predictions(predictions_dir: Path | None, predictions_glob: str) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    if predictions_dir is not None:
        paths = sorted(predictions_dir.rglob("porter6_predictions_merged.csv"))
        if not paths:
            raise SystemExit(f"No porter6_predictions_merged.csv under {predictions_dir}")
        for path in paths:
            frames.append(pd.read_csv(path))
    elif predictions_glob:
        paths = sorted(glob(predictions_glob))
        if not paths:
            raise SystemExit(f"No files matched: {predictions_glob}")
        for path in paths:
            frames.append(pd.read_csv(path))
    else:
        raise SystemExit("Provide --predictions-glob or --predictions-dir")

    preds = pd.concat(frames, ignore_index=True)
    if "id" not in preds.columns:
        raise SystemExit("Predictions must include an 'id' column")
    preds["id"] = preds["id"].astype(str).str.strip()
    preds = preds.drop_duplicates(subset=["id"], keep="first")
    return preds


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--map-tsv", type=Path, required=True, help="orig_to_uniq.tsv from pangenome_dedupe.py")
    ap.add_argument(
        "--predictions-glob",
        type=str,
        default="",
        help="Glob for porter6_predictions_merged.csv (e.g. 'results/batch_*/porter6_predictions_merged.csv')",
    )
    ap.add_argument(
        "--predictions-dir",
        type=Path,
        default=None,
        help="If set, load every .../porter6_predictions_merged.csv under this directory",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=Path("recombined"),
        help="Output directory for per-genome files (default: recombined)",
    )
    args = ap.parse_args()

    mapping = pd.read_csv(args.map_tsv, sep="\t")
    required = {"uniq_id", "header_id", "source_file"}
    if not required.issubset(mapping.columns):
        raise SystemExit(f"Unexpected columns in map: {list(mapping.columns)} (need {sorted(required)})")

    preds = load_predictions(args.predictions_dir, args.predictions_glob)
    merged = mapping.merge(preds, left_on="uniq_id", right_on="id", how="left")
    merged["genome_id"] = merged["source_file"].map(genome_id_from_source)

    if "q3" not in merged.columns or "q8" not in merged.columns:
        raise SystemExit("Predictions must contain q3 and q8 columns")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    missing = int(merged["q3"].isna().sum())
    if missing:
        print(f"Warning: {missing} rows missing predictions (check uniq_id / batch runs)")

    for genome_id, grp in merged.groupby("genome_id", sort=True):
        q3_out = args.out_dir / f"{genome_id}.q3.csv"
        q8_out = args.out_dir / f"{genome_id}.q8.csv"

        q3_df = grp[["header_id", "q3"]].rename(columns={"header_id": "id", "q3": "q3_decoded"})
        q8_df = grp[["header_id", "q8"]].rename(columns={"header_id": "id", "q8": "q8_decoded"})
        q3_df.to_csv(q3_out, index=False)
        q8_df.to_csv(q8_out, index=False)
        print(f"Wrote {q3_out} ({len(q3_df)} rows)")
        print(f"Wrote {q8_out} ({len(q8_df)} rows)")


if __name__ == "__main__":
    main()
