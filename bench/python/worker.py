#!/usr/bin/env python
"""Concurrency worker (Python): one simulated session issuing a mixed query stream.
  python -I bench/python/worker.py --scale 1M --impl raw_anndata_py --worker 1 --duration 60 --go <barrier> [--run-id ID] [--sessions N]
"""
from __future__ import annotations

import argparse
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run as H  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scale", required=True); ap.add_argument("--impl", required=True)
    ap.add_argument("--worker", type=int, default=1); ap.add_argument("--duration", type=float, default=60)
    ap.add_argument("--warmup", type=float, default=10); ap.add_argument("--go", default=""); ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--sessions", dest="n", type=int, default=1); ap.add_argument("--run-id", default="")
    ap.add_argument("--config", default=os.path.join(H.BENCH, "config", "bench.toml"))
    a = ap.parse_args()
    opts = argparse.Namespace(scale=a.scale, config=a.config, cache="warm", only=a.impl, resume=False, run_id=a.run_id, single="",
                              rep=1, no_tracemalloc=True, start_time=None)
    ctx = H.Ctx(opts)
    ctx.open_output(f"python_conc_{a.impl}_{ctx.scale}_process_n{a.n}_w{a.worker}")
    obj = H.IMPLS[a.impl](ctx)
    all_ids = obj.cells
    rng = random.Random(a.seed + a.worker)
    mix = ctx.cfg["concurrency"]["mix"]; names = list(mix); weights = [mix[k] for k in names]

    def pick():
        name = rng.choices(names, weights)[0]
        if name == "feature1_all":
            return dict(workload="feature", layer="data", gene_set=f"1:{rng.randint(0, 2)}", n_features=1, cell_subset="all", extra="")
        if name == "feature10_rand10":
            return dict(workload="feature", layer="data", gene_set=f"10:{rng.randint(0, 2)}", n_features=10, cell_subset="rand10", extra="")
        if name == "meta_col":
            return dict(workload="metadata", layer="", gene_set="", n_features=0, cell_subset="all", extra=rng.choice(["cell_type", "nCount_RNA"]))
        if name == "reduction":
            return dict(workload="reduction", layer="", gene_set="", n_features=0, cell_subset="all", extra="X_umap|1-2")
        return dict(workload="fetch_cells", layer="", gene_set="", n_features=0, cell_subset="", extra=rng.choice(ctx.cfg["workloads"]["fetch_cells"]))

    t_warm = time.time()
    while (not os.path.exists(a.go)) if a.go else (time.time() - t_warm < a.warmup):
        w = pick()
        if H.run_work(ctx, obj, w, all_ids) is None:
            continue
        if time.time() - t_warm > 300:
            break
    t_end = time.time() + a.duration; n = 0
    while time.time() < t_end:
        w = pick()
        res, wall, cpu, alloc = H.timed(lambda: H.run_work(ctx, obj, w, all_ids), False)
        if res is None:
            continue
        chk, nout, exp = res
        ok = H.checksum_ok(chk, exp) if rng.random() < ctx.cfg["concurrency"]["checksum_fraction"] else ""
        ctx.row(implementation=a.impl, storage=H.STORAGE[a.impl], workload=H.label(w), n_features=w["n_features"], gene_set=H.gene_set_index(w),
                cell_subset=w["cell_subset"], layer=w["layer"], n_cells_out=nout, phase="concurrent", rep=a.worker, wall_ms=wall, cpu_ms=cpu,
                result_checksum=repr(chk) if isinstance(chk, float) else str(chk), checksum_ok=ok, note=f"mode=process;n={a.n}", threads=a.n)
        n += 1
    ctx.close()
    print(f"worker {a.worker}: {n} queries")


if __name__ == "__main__":
    main()
