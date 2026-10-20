#!/usr/bin/env python
"""Cold-cache runs: evict, then run ONE workload in a fresh process, repeat.

  python -I tools/run_cold.py --scale 1M --lang julia --impl jl_h5ad_csc [--workloads "feature|data|1:0|all,metadata|full"]
      [--reps 2] [--run-id ID] [--config config/bench.toml] [--engine r]

Each repetition: tools/cold.sh evicts the storage's files, then the language harness is started with
`--cache cold --single <workload> --rep <k>` so it opens the source and runs the workload once. If
eviction could not be verified the rows are written with cache_state=cold_unverified.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
import tomllib

BENCH = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
STORAGE = {"jl_h5ad_csc": "h5ad_csc", "jl_h5ad_csr": "h5ad_csr", "jl_bpcells": "bpcells", "r_seurat_bpcells": "seurat",
           "r_sce_backed": "sce", "raw_bpcells_r": "bpcells", "raw_anndata_py": "h5ad_csc", "raw_h5py_py": "h5ad_csc",
           "py_scuba_sc": "h5ad_csc"}
DEFAULT_WORKLOADS = ["feature|data|1:0|all", "feature|data|10:0|all", "feature|data|100:0|all", "feature|data|10:0|rand10",
                     "metadata|full", "reduction|X_umap|1-2"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scale", required=True)
    ap.add_argument("--lang", required=True, choices=["julia", "r", "python"])
    ap.add_argument("--impl", required=True)
    ap.add_argument("--engine", default="r")
    ap.add_argument("--workloads", default=",".join(DEFAULT_WORKLOADS))
    ap.add_argument("--reps", type=int, default=None)
    ap.add_argument("--run-id", default="")
    ap.add_argument("--config", default=os.path.join(BENCH, "config", "bench.toml"))
    args = ap.parse_args()
    with open(args.config, "rb") as fh:
        cfg = tomllib.load(fh)
    reps = args.reps or cfg["reps"]["cold"].get(args.scale, 1)
    data = os.path.join(os.path.expanduser(cfg["suite"]["data_dir"]), args.scale)
    run_id = args.run_id or cfg["suite"]["run_id"] or time.strftime("%Y%m%dT%H%M%S", time.gmtime())
    timeout = cfg["reps"]["timeout_s"]
    for w in args.workloads.split(","):
        for rep in range(1, reps + 1):
            ev = subprocess.run(["bash", os.path.join(BENCH, "tools", "cold.sh"), data, STORAGE[args.impl]], capture_output=True, text=True)
            cache = "cold" if ev.returncode == 0 else "cold_unverified"
            print(f"[cold] {args.impl} {w} rep {rep}: {ev.stdout.strip().splitlines()[-1] if ev.stdout.strip() else ev.stderr.strip()} -> {cache}", flush=True)
            common = ["--scale", args.scale, "--run-id", run_id, "--cache", cache, "--single", w, "--rep", str(rep), "--only", args.impl,
                      "--start-time", f"{time.time():.3f}"]
            if args.lang == "julia":
                cmd = ["julia", f"--project={BENCH}", "-t", "1", os.path.join(BENCH, "julia", "run.jl"), *common]
            elif args.lang == "r":
                cmd = ["Rscript", os.path.join(BENCH, "r", "run.R"), *common, "--engine", args.engine]
            else:
                cmd = ["conda", "run", "-n", cfg["suite"]["conda_env"], "--no-capture-output", "python", "-I",
                       os.path.join(BENCH, "python", "run.py"), *common]
            try:
                r = subprocess.run(cmd, timeout=timeout)
                if r.returncode != 0:
                    print(f"  harness exited {r.returncode}", file=sys.stderr)
            except subprocess.TimeoutExpired:
                print(f"  TIMEOUT after {timeout}s", file=sys.stderr)


if __name__ == "__main__":
    main()
