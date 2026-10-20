#!/usr/bin/env python
"""Pre-run checks and environment capture for the SCUBA benchmark suite.

  python -I tools/preflight.py --scale 1M [--config config/bench.toml] [--run-id ID] [--no-env]

Checks free disk and RAM against the projected dataset size, the presence of julia / Rscript /
vmtouch / the conda env, and required R packages; then writes results/env/<run_id>.json with
tool and package versions, git SHAs and hardware facts. Exit code 1 on a hard failure.
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import time
import tomllib

BENCH = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def sh(cmd, **kw):
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=kw.get("timeout", 600), check=False).stdout.strip()
    except Exception as e:  # noqa: BLE001
        return f"<error: {e}>"


def git(path, *args):
    return sh(["git", "-C", os.path.expanduser(path), *args])


def projected_bytes(sc):
    nnz = sc["cells"] * sc["nnz_per_cell"]
    per_layer = nnz * 8                               # float32 data + int32 indices
    csr_layers = len(sc.get("csr_layers", ["X", "counts"]))
    h5ad = per_layer * (2 + csr_layers)
    bpcells = nnz * 7                                 # packed uint counts (~3 B) + packed double data (~4 B)
    return {"nnz": nnz, "one_layer": per_layer, "h5ad": h5ad, "bpcells": bpcells, "total": h5ad + bpcells,
            "ram_generation": int(per_layer * 3.5)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scale", default=None)
    ap.add_argument("--config", default=os.path.join(BENCH, "config", "bench.toml"))
    ap.add_argument("--run-id", default=None)
    ap.add_argument("--no-env", action="store_true")
    args = ap.parse_args()
    with open(args.config, "rb") as fh:
        cfg = tomllib.load(fh)
    data_root = os.path.expanduser(cfg["suite"]["data_dir"])
    os.makedirs(data_root, exist_ok=True)
    failures, warnings = [], []

    # disk and RAM
    free = shutil.disk_usage(data_root).free
    ram = int(sh(["sysctl", "-n", "hw.memsize"]) or 0) if platform.system() == "Darwin" else os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    print(f"free disk at {data_root}: {free / 2**30:.0f} GiB; RAM: {ram / 2**30:.0f} GiB")
    if args.scale:
        sc = cfg["scales"][args.scale]
        p = projected_bytes(sc)
        done = os.path.exists(os.path.join(data_root, args.scale, "MANIFEST.json"))
        print(f"{args.scale}: projected {p['total'] / 2**30:.0f} GiB on disk, {p['ram_generation'] / 2**30:.0f} GiB RAM to generate"
              + (" (dataset already present)" if done else ""))
        if not done and free < p["total"] * 1.2:
            failures.append(f"insufficient disk for {args.scale}: need {p['total'] * 1.2 / 2**30:.0f} GiB")
        if not done and ram < p["ram_generation"]:
            failures.append(f"insufficient RAM to generate {args.scale}: need {p['ram_generation'] / 2**30:.0f} GiB")

    # tools
    for tool in ("julia", "Rscript", "conda"):
        if shutil.which(tool) is None:
            failures.append(f"{tool} not on PATH")
    if shutil.which("vmtouch") is None:
        warnings.append("vmtouch not found (brew install vmtouch); cold-cache runs will be marked cold_unverified")
    env = cfg["suite"]["conda_env"]
    envs = sh(["conda", "env", "list"])
    if env not in envs:
        failures.append(f"conda env {env} not found")
    rpk = sh(["Rscript", "-e", 'for (p in c("SeuratObject","Seurat","BPCells","HDF5Array","SingleCellExperiment","JuliaCall","arrow","ps","RcppTOML","jsonlite","SCUBA")) cat(p, "=", if (requireNamespace(p, quietly=TRUE)) as.character(packageVersion(p)) else "MISSING", "\\n")'])
    missing_r = [l.split(" = ")[0] for l in rpk.splitlines() if l.endswith("MISSING")]
    if missing_r:
        failures.append("missing R packages: " + ", ".join(missing_r))
    anndata_r = sh(["Rscript", "-e", 'cat(if (requireNamespace("anndata", quietly=TRUE)) as.character(packageVersion("anndata")) else "MISSING")'])
    if anndata_r != "MISSING" and tuple(int(x) for x in anndata_r.split(".")[:2]) >= (0, 8):
        warnings.append(f"R anndata {anndata_r} returns python-class objects; the legacy r_anndata_backed path will be skipped "
                        "(remotes::install_version('anndata', '0.7.5.6') to enable it)")

    for w in warnings:
        print("WARNING:", w)
    for f in failures:
        print("FAIL:", f)

    if not args.no_env:
        run_id = args.run_id or cfg["suite"].get("run_id") or time.strftime("%Y%m%dT%H%M%S")
        env_dir = os.path.join(BENCH, cfg["suite"]["results_dir"], "env")
        os.makedirs(env_dir, exist_ok=True)
        info = {
            "run_id": run_id, "suite_version": cfg["suite"]["version"], "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "host": platform.node(), "platform": platform.platform(),
            "hardware": {k: sh(["sysctl", "-n", k]) for k in ("hw.model", "hw.memsize", "hw.ncpu", "machdep.cpu.brand_string")} if platform.system() == "Darwin" else {},
            "os": sh(["sw_vers"]).splitlines() if platform.system() == "Darwin" else platform.version(),
            "julia": {"version": sh(["julia", "--version"]),
                      "status": sh(["julia", f"--project={BENCH}", "-e", "using Pkg; Pkg.status(io=stdout)"])},
            "scuba_jl": {"sha": git(os.path.join(BENCH, ".."), "rev-parse", "HEAD"), "dirty": bool(git(os.path.join(BENCH, ".."), "status", "--porcelain"))},
            "scuba_r": {"path": cfg["suite"]["r_scuba_path"], "sha": git(cfg["suite"]["r_scuba_path"], "rev-parse", "HEAD"),
                        "branch": git(cfg["suite"]["r_scuba_path"], "rev-parse", "--abbrev-ref", "HEAD"),
                        "dirty": bool(git(cfg["suite"]["r_scuba_path"], "status", "--porcelain"))},
            "R": {"version": sh(["Rscript", "-e", "cat(R.version.string)"]), "packages": rpk.splitlines(), "anndata": anndata_r},
            "python": {"env": env, "packages": [p for p in sh(["conda", "list", "-n", env, "--json"]).split("\n")][:1]},
            "tools": {"vmtouch": shutil.which("vmtouch")},
            "warnings": warnings, "failures": failures,
        }
        try:
            pk = json.loads(sh(["conda", "list", "-n", env, "--json"]))
            info["python"]["packages"] = {p["name"]: p["version"] for p in pk if p["name"] in
                                          ("python", "anndata", "h5py", "scipy", "numpy", "pandas", "pyarrow", "psutil", "juliacall", "zarr")}
        except Exception:
            pass
        path = os.path.join(env_dir, f"{run_id}.json")
        with open(path, "w") as fh:
            json.dump(info, fh, indent=2)
        print(f"environment written to {path}")
        print(run_id)
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
