#!/usr/bin/env python
"""Merge raw per-process CSVs into results/results.parquet and results/results.csv.

  python -I tools/merge_results.py [--results-dir bench/results] [--run-id ID]
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import pandas as pd

BENCH = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
NUMERIC = ["n_features", "gene_set", "n_cells_out", "rep", "threads", "wall_ms", "cpu_ms", "open_ms", "rss_mb",
           "maxrss_mb", "alloc_mb", "gc_ms"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", default=os.path.join(BENCH, "results"))
    ap.add_argument("--run-id", default=None)
    args = ap.parse_args()
    pattern = os.path.join(args.results_dir, "raw", args.run_id or "*", "*.csv")
    files = sorted(glob.glob(pattern))
    if not files:
        print("no raw result files found at", pattern)
        sys.exit(1)
    frames = []
    for f in files:
        try:
            df = pd.read_csv(f, dtype=str, keep_default_na=False)
        except pd.errors.EmptyDataError:
            continue
        df["source_file"] = os.path.relpath(f, args.results_dir)
        frames.append(df)
    all_df = pd.concat(frames, ignore_index=True)
    for c in NUMERIC:
        if c in all_df:
            all_df[c] = pd.to_numeric(all_df[c], errors="coerce")
    if "checksum_ok" in all_df:
        all_df["checksum_ok"] = all_df["checksum_ok"].map({"true": True, "false": False, "TRUE": True, "FALSE": False, "True": True, "False": False}).astype("boolean")
    out_parquet = os.path.join(args.results_dir, "results.parquet")
    out_csv = os.path.join(args.results_dir, "results.csv")
    all_df.to_parquet(out_parquet, index=False)
    all_df.to_csv(out_csv, index=False)
    bad = all_df[(all_df["status"] == "ok") & (all_df["checksum_ok"] == False)]  # noqa: E712
    print(f"merged {len(files)} files, {len(all_df)} rows -> {out_parquet}")
    print(all_df.groupby(["language", "implementation", "scale", "cache_state"]).size().to_string())
    if len(bad):
        print(f"WARNING: {len(bad)} rows with checksum mismatches")
        print(bad[["implementation", "scale", "workload", "n_features", "gene_set", "cell_subset", "result_checksum"]].head(20).to_string())


if __name__ == "__main__":
    main()
