#!/usr/bin/env python
"""Concurrency driver: N simultaneous sessions against one dataset, multi-process (any language) or Julia threads.

  python -I tools/concurrency.py --scale 1M --lang julia --impl jl_h5ad_csc --sessions 1,2,4,8,16 [--run-id ID]
  python -I tools/concurrency.py --scale 1M --lang julia --impl jl_bpcells --mode threads --share shared --sessions 1,4,16
  python -I tools/concurrency.py --scale 1M --lang r --impl r_seurat_bpcells --engine julia --sessions 1,4
  python -I tools/concurrency.py --scale 1M --lang python --impl raw_anndata_py --sessions 1,4

Workers warm up until a barrier file appears (created once all workers have started plus
[concurrency].warmup_s), then run for [concurrency].duration_s. Aggregates per N into
results/raw/<run_id>/concurrency_summary.csv (throughput, p50/p95/p99 per workload).
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import statistics
import subprocess
import tempfile
import time
import tomllib

BENCH = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def worker_cmd(lang, cfg, common, engine):
    if lang == "julia":
        return ["julia", f"--project={BENCH}", "-t", "1", os.path.join(BENCH, "julia", "worker.jl"), *common]
    if lang == "r":
        return ["Rscript", os.path.join(BENCH, "r", "worker.R"), *common, "--engine", engine]
    return ["conda", "run", "-n", cfg["suite"]["conda_env"], "--no-capture-output", "python", "-I",
            os.path.join(BENCH, "python", "worker.py"), *common]


def summarize(run_dir, pattern, n, lang, impl, mode, share, duration, out_rows):
    rows = []
    for f in glob.glob(os.path.join(run_dir, pattern)):
        with open(f) as fh:
            rows += [r for r in csv.DictReader(fh) if r["phase"] == "concurrent"]
    if not rows:
        return
    walls = [float(r["wall_ms"]) for r in rows if r["wall_ms"]]
    def pct(v, q):
        v = sorted(v); return v[min(len(v) - 1, int(q * len(v)))]
    by_w = {}
    for r in rows:
        if r["wall_ms"]:
            by_w.setdefault(r["workload"], []).append(float(r["wall_ms"]))
    bad = sum(1 for r in rows if r["checksum_ok"] == "false")
    out_rows.append(dict(language=lang, implementation=impl, mode=mode, share=share, n=n, queries=len(rows), duration_s=duration,
                         throughput_qps=len(rows) / duration, p50_ms=statistics.median(walls), p95_ms=pct(walls, 0.95), p99_ms=pct(walls, 0.99),
                         workload="all", checksum_mismatches=bad))
    for w, v in sorted(by_w.items()):
        out_rows.append(dict(language=lang, implementation=impl, mode=mode, share=share, n=n, queries=len(v), duration_s=duration,
                             throughput_qps=len(v) / duration, p50_ms=statistics.median(v), p95_ms=pct(v, 0.95), p99_ms=pct(v, 0.99),
                             workload=w, checksum_mismatches=0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scale", required=True); ap.add_argument("--lang", required=True, choices=["julia", "r", "python"])
    ap.add_argument("--impl", required=True); ap.add_argument("--engine", default="r")
    ap.add_argument("--mode", default="process", choices=["process", "threads"]); ap.add_argument("--share", default="shared")
    ap.add_argument("--sessions", default=None); ap.add_argument("--run-id", default="")
    ap.add_argument("--config", default=os.path.join(BENCH, "config", "bench.toml"))
    a = ap.parse_args()
    with open(a.config, "rb") as fh:
        cfg = tomllib.load(fh)
    conc = cfg["concurrency"]
    ns = [int(x) for x in (a.sessions.split(",") if a.sessions else conc["n"])]
    run_id = a.run_id or cfg["suite"]["run_id"] or time.strftime("%Y%m%dT%H%M%S", time.gmtime())
    run_dir = os.path.join(BENCH, cfg["suite"]["results_dir"], "raw", run_id)
    os.makedirs(run_dir, exist_ok=True)
    out_rows = []
    for n in ns:
        # remove stale worker files for this configuration
        for f in glob.glob(os.path.join(run_dir, f"{a.lang}_conc_{a.impl}_*{a.scale}_{a.mode}_n{n}_*.csv")):
            os.remove(f)
        go = os.path.join(tempfile.mkdtemp(prefix="scuba-go-"), "go")
        common = ["--scale", a.scale, "--impl", a.impl, "--duration", str(conc["duration_s"]), "--warmup", str(conc["warmup_s"]),
                  "--go", go, "--run-id", run_id, "--sessions", str(n)]
        procs = []
        t0 = time.time()
        if a.mode == "threads":
            cmd = ["julia", f"--project={BENCH}", "-t", str(n), os.path.join(BENCH, "julia", "worker.jl"), *common, "--mode", "threads", "--share", a.share]
            procs.append(subprocess.Popen(cmd))
        else:
            for w in range(1, n + 1):
                procs.append(subprocess.Popen(worker_cmd(a.lang, cfg, common + ["--worker", str(w)], a.engine)))
        print(f"[conc] {a.lang} {a.impl} n={n} mode={a.mode}: {len(procs)} process(es) started; warm-up {conc['warmup_s']}s", flush=True)
        time.sleep(conc["warmup_s"] + (60 if a.lang != "python" else 20))   # generous: JIT / readRDS / JuliaCall init before the barrier
        open(go, "w").close()
        for p in procs:
            p.wait()
        print(f"[conc] n={n} done in {time.time() - t0:.0f}s", flush=True)
        summarize(run_dir, f"{a.lang}_conc_{a.impl}_{a.scale}_{a.mode}_n{n}_*.csv", n, a.lang, a.impl, a.mode, a.share if a.mode == "threads" else "process",
                  conc["duration_s"], out_rows)
    path = os.path.join(run_dir, "concurrency_summary.csv")
    exists = os.path.exists(path)
    with open(path, "a", newline="") as fh:
        wri = csv.DictWriter(fh, fieldnames=list(out_rows[0].keys()) if out_rows else ["language"])
        if not exists:
            wri.writeheader()
        for r in out_rows:
            wri.writerow(r)
    for r in out_rows:
        if r["workload"] == "all":
            print(f"  n={r['n']}: {r['throughput_qps']:.1f} q/s, p50 {r['p50_ms']:.1f} ms, p95 {r['p95_ms']:.1f} ms")


if __name__ == "__main__":
    main()
