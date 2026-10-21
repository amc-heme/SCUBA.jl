#!/usr/bin/env python
"""Python benchmark harness: raw anndata / h5py baselines and the scuba-sc package.

  conda run -n scuba-fixtures python -I bench/python/run.py --scale 10k [--cache warm|cold] [--only raw_anndata_py,...]
      [--resume] [--run-id ID] [--single "feature|data|10:0|all" --rep N] [--no-tracemalloc] [--start-time SECS]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import resource
import socket
import sys
import time
import tomllib
import tracemalloc

T_PROCESS_START = time.time()
import numpy as np  # noqa: E402

BENCH = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
SCHEMA = ["suite_version", "run_id", "timestamp", "host", "language", "implementation", "engine", "storage", "scale",
          "workload", "n_features", "gene_set", "cell_subset", "layer", "n_cells_out", "cache_state", "phase", "rep",
          "threads", "wall_ms", "cpu_ms", "open_ms", "rss_mb", "maxrss_mb", "alloc_mb", "gc_ms", "result_checksum",
          "checksum_ok", "status", "note"]
STORAGE = {"raw_anndata_py": "h5ad_csc", "raw_h5py_py": "h5ad_csc", "py_scuba_sc": "h5ad_csc"}


def sha1_lines(v):
    return hashlib.sha1("\n".join(map(str, v)).encode()).hexdigest()


def rss_mb():
    try:
        import psutil
        return psutil.Process().memory_info().rss / 2**20
    except Exception:
        return float("nan")


def maxrss_mb():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / 2**20 if platform.system() == "Darwin" else r / 1024


def checksum_ok(got, expected, rtol=1e-5):
    if expected is None:
        return ""
    if isinstance(got, (int, float)) and isinstance(expected, (int, float)):
        return str(abs(got - expected) <= rtol * abs(expected) + 1e-9).lower()
    a, b = str(got).split(":", 1), str(expected).split(":", 1)
    if len(a) == 2 and len(b) == 2:
        if a[0] != b[0]:
            return "false"
        try:
            return str(abs(float(a[1]) - float(b[1])) <= rtol * abs(float(b[1])) + 1e-9).lower()
        except ValueError:
            return str(a[1] == b[1]).lower()
    return str(got == expected).lower()


class Ctx:
    def __init__(self, opts):
        with open(opts.config, "rb") as fh:
            self.cfg = tomllib.load(fh)
        self.scale = opts.scale
        self.data = os.path.join(os.path.expanduser(self.cfg["suite"]["data_dir"]), opts.scale)
        with open(os.path.join(self.data, "MANIFEST.json")) as fh:
            self.manifest = json.load(fh)
        assert self.manifest.get("complete"), "dataset incomplete"
        with open(os.path.join(self.data, "selection.json")) as fh:
            self.selection = json.load(fh)
        with open(os.path.join(self.data, "expected.json")) as fh:
            self.expected = json.load(fh)
        self.run_id = opts.run_id or self.cfg["suite"]["run_id"] or time.strftime("%Y%m%dT%H%M%S", time.gmtime())
        self.out_dir = os.path.join(BENCH, self.cfg["suite"]["results_dir"], "raw", self.run_id)
        os.makedirs(self.out_dir, exist_ok=True)
        self.cache_state = opts.cache
        self.host = socket.gethostname()
        self.fh = None
        self.existing = set()
        self.opts = opts

    def open_output(self, name):
        path = os.path.join(self.out_dir, name + ".csv")
        fresh = not os.path.exists(path) or os.path.getsize(path) == 0
        self.existing = set() if fresh else existing_keys(path)
        self.fh = open(path, "a")
        if fresh:
            self.fh.write(",".join(SCHEMA) + "\n")
            self.fh.flush()

    def close(self):
        if self.fh:
            self.fh.close(); self.fh = None

    def row(self, **kw):
        vals = {"suite_version": self.cfg["suite"]["version"], "run_id": self.run_id,
                "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "host": self.host, "language": "python",
                "engine": "none", "scale": self.scale, "n_features": 0, "gene_set": -1, "cell_subset": "", "layer": "",
                "n_cells_out": 0, "cache_state": self.cache_state, "phase": "fetch", "rep": 0, "threads": 1,
                "wall_ms": "", "cpu_ms": "", "open_ms": 0, "rss_mb": round(rss_mb(), 2), "maxrss_mb": round(maxrss_mb(), 2),
                "alloc_mb": "", "gc_ms": "", "result_checksum": "", "checksum_ok": "", "status": "ok", "note": ""}
        vals.update(kw)
        out = []
        for c in SCHEMA:
            s = vals.get(c, "")
            s = "" if s is None else (f"{s:.4f}" if isinstance(s, float) else str(s))
            out.append('"' + s.replace('"', '""') + '"' if any(ch in s for ch in ',"\n') else s)
        self.fh.write(",".join(out) + "\n"); self.fh.flush()
        self.existing.add(row_key(vals["implementation"], vals["workload"], vals["n_features"], vals["gene_set"],
                                  vals["cell_subset"], vals["layer"], self.cache_state, vals["phase"], vals["rep"]))


def row_key(*parts):
    return "|".join(str(p) for p in parts)


def existing_keys(path):
    import csv
    keys = set()
    with open(path) as fh:
        for r in csv.DictReader(fh):
            keys.add(row_key(r["implementation"], r["workload"], r["n_features"], r["gene_set"], r["cell_subset"],
                             r["layer"], r["cache_state"], r["phase"], r["rep"]))
    return keys


def timed(f, use_tracemalloc):
    if use_tracemalloc:
        tracemalloc.start()
    c0 = time.process_time_ns(); t0 = time.perf_counter_ns()
    v = f()
    t1 = time.perf_counter_ns(); c1 = time.process_time_ns()
    alloc = ""
    if use_tracemalloc:
        _, peak = tracemalloc.get_traced_memory(); tracemalloc.stop(); alloc = peak / 2**20
    return v, (t1 - t0) / 1e6, (c1 - c0) / 1e6, alloc


# ---------------------------------------------------------------------------
# implementations

class RawAnnData:
    def __init__(self, ctx):
        import anndata as ad
        self.adata = ad.read_h5ad(os.path.join(ctx.data, "csc.h5ad"), backed="r")
        self.cells = np.asarray(self.adata.obs_names.astype(str))
        self.pos = {c: i for i, c in enumerate(self.cells)}

    def feature(self, genes, cells, layer):
        import scipy.sparse as sp
        sub = self.adata[:, genes] if cells is None else self.adata[[self.pos[c] for c in cells], genes]
        X = sub.X if layer == "data" else sub.layers["counts"]
        if hasattr(X, "to_memory"):
            X = X.to_memory()
        if sp.issparse(X):
            return float(X.sum(dtype=np.float64)), X.shape[0]
        return float(np.asarray(X, dtype=np.float64).sum()), X.shape[0]

    def metadata(self, col):
        obs = self.adata.obs
        if col == "full":
            return f"{obs.shape[0]}:{float(obs['nCount_RNA'].sum())!r}", obs.shape[0]
        v = obs[col]
        if v.dtype.kind in "fiu":
            return f"{len(v)}:{float(v.sum())!r}", len(v)
        return f"{len(v)}:{sha1_lines(v.astype(str).tolist())}", len(v)

    def reduction(self, name, lo, hi):
        arr = np.asarray(self.adata.obsm[name])[:, lo - 1:hi]
        return float(arr.astype(np.float64).sum()), arr.shape[0]

    def fetch_cells(self, col, level):
        ids = sorted(self.cells[(self.adata.obs[col].astype(str) == level).to_numpy()].tolist())
        return f"{len(ids)}:{sha1_lines(ids)}", len(ids)

    def unique_values(self, col):
        lv = sorted(set(self.adata.obs[col].astype(str).tolist()))
        return f"{len(lv)}:{sha1_lines(lv)}", len(lv)


class RawH5py:
    """Direct CSC slicing with h5py: the floor for the HDF5 path (features only)."""
    def __init__(self, ctx):
        import h5py
        self.f = h5py.File(os.path.join(ctx.data, "csc.h5ad"), "r")
        self.var = np.asarray(self.f["var"]["gene"][:]).astype(str)
        self.gpos = {g: i for i, g in enumerate(self.var)}
        self.cells = np.asarray(self.f["obs"]["cell_id"][:]).astype(str)
        self.cpos = {c: i for i, c in enumerate(self.cells)}
        self.indptr = {"data": self.f["X"]["indptr"][:], "counts": self.f["layers"]["counts"]["indptr"][:]}
        self.n_obs = len(self.cells)

    def feature(self, genes, cells, layer):
        g = self.f["X"] if layer == "data" else self.f["layers"]["counts"]
        indptr = self.indptr[layer]
        mask = None
        if cells is not None:
            mask = np.zeros(self.n_obs, dtype=bool); mask[[self.cpos[c] for c in cells]] = True
        total = 0.0
        for gene in genes:
            j = self.gpos[gene]
            a, b = int(indptr[j]), int(indptr[j + 1])
            if b <= a:
                continue
            d = g["data"][a:b]
            if mask is None:
                total += float(d.astype(np.float64).sum())
            else:
                idx = g["indices"][a:b]
                total += float(d[mask[idx]].astype(np.float64).sum())
        return total, self.n_obs if cells is None else len(cells)


class ScubaSC:
    def __init__(self, ctx):
        import scuba
        self.scuba = scuba
        self.ds = scuba.open(os.path.join(ctx.data, "csc.h5ad"))
        self.cells = scuba.get_all_cells(self.ds)

    def feature(self, genes, cells, layer):
        df = self.scuba.fetch_data(self.ds, [f"X_{g}" for g in genes], layer=None if layer == "data" else layer, cells=cells)
        return float(df.to_numpy(dtype=np.float64).sum()), df.shape[0]

    def metadata(self, col):
        if col == "full":
            df = self.scuba.fetch_metadata(self.ds)
            return f"{df.shape[0]}:{float(df['nCount_RNA'].sum())!r}", df.shape[0]
        v = self.scuba.fetch_metadata(self.ds, col, return_class="vector")
        if str(v.dtype).startswith(("float", "int", "Int", "Float")):
            return f"{len(v)}:{float(v.astype(float).sum())!r}", len(v)
        return f"{len(v)}:{sha1_lines(v.astype(str).tolist())}", len(v)

    def reduction(self, name, lo, hi):
        df = self.scuba.fetch_reduction(self.ds, name, dims=list(range(lo, hi + 1)))
        return float(df.to_numpy(dtype=np.float64).sum()), df.shape[0]

    def fetch_cells(self, col, level):
        ids = sorted(self.scuba.fetch_cells(self.ds, col, [level]))
        return f"{len(ids)}:{sha1_lines(ids)}", len(ids)

    def unique_values(self, col):
        lv = sorted(self.scuba.unique_values(self.ds, col))
        return f"{len(lv)}:{sha1_lines(lv)}", len(lv)


IMPLS = {"raw_anndata_py": RawAnnData, "raw_h5py_py": RawH5py, "py_scuba_sc": ScubaSC}


# ---------------------------------------------------------------------------
# workloads (same enumeration as the Julia harness)

def workloads(ctx):
    wl = ctx.cfg["workloads"]
    W = []
    keys = sorted(ctx.selection["gene_sets"], key=lambda k: tuple(int(x) for x in k.split(":")))
    for key in keys:
        n, rep = (int(x) for x in key.split(":"))
        layers = list(wl["layers_all_sets"]) + (list(wl["layers_first_set"]) if rep == 0 else [])
        for layer in layers:
            for sub in wl["cell_subsets"]:
                W.append(dict(workload="feature", layer=layer, gene_set=key, n_features=n, cell_subset=sub, extra=""))
    for m in wl["metadata"]:
        W.append(dict(workload="metadata", layer="", gene_set="", n_features=0, cell_subset="all", extra=m))
    for r in wl["reductions"]:
        W.append(dict(workload="reduction", layer="", gene_set="", n_features=0, cell_subset="all", extra=f"{r['name']}|{r['dims'][0]}-{r['dims'][1]}"))
    for f in wl["fetch_cells"]:
        W.append(dict(workload="fetch_cells", layer="", gene_set="", n_features=0, cell_subset="", extra=f))
    for u in wl["unique_values"]:
        W.append(dict(workload="unique_values", layer="", gene_set="", n_features=0, cell_subset="", extra=u))
    return W


def parse_single(spec):
    p = spec.split("|")
    if p[0] in ("feature", "matrix_sparse"):
        return dict(workload=p[0], layer=p[1], gene_set=p[2], n_features=int(p[2].split(":")[0]), cell_subset=p[3], extra="")
    if p[0] in ("metadata", "reduction"):
        return dict(workload=p[0], layer="", gene_set="", n_features=0, cell_subset="all", extra="|".join(p[1:]))
    return dict(workload=p[0], layer="", gene_set="", n_features=0, cell_subset="", extra=p[1])


def cells_for_subset(ctx, all_ids, sub):
    if sub == "all":
        return None
    s = ctx.selection["cell_subsets"][sub]
    if sub.startswith("contig"):
        return [str(c) for c in all_ids[s[0]:s[1]]]
    return [str(all_ids[i]) for i in s]


def run_work(ctx, impl_obj, w, all_ids):
    k = w["workload"]
    if k == "feature":
        genes = ctx.selection["gene_sets"][w["gene_set"]]
        cells = cells_for_subset(ctx, all_ids, w["cell_subset"])
        chk, n = impl_obj.feature(genes, cells, w["layer"])
        e = ctx.expected["feature"].get(f"{w['layer']}|{w['gene_set']}|{w['cell_subset']}")
        return chk, n, (e or {}).get("sum")
    if not hasattr(impl_obj, k):
        return None
    if k == "metadata":
        chk, n = impl_obj.metadata(w["extra"]); return chk, n, ctx.expected["metadata"][w["extra"]]
    if k == "reduction":
        name, dims = w["extra"].split("|"); lo, hi = (int(x) for x in dims.split("-"))
        chk, n = impl_obj.reduction(name, lo, hi); return chk, n, ctx.expected["reduction"][w["extra"]]
    if k == "fetch_cells":
        col, level = w["extra"].split(":"); chk, n = impl_obj.fetch_cells(col, level); return chk, n, ctx.expected["fetch_cells"][w["extra"]]
    if k == "unique_values":
        chk, n = impl_obj.unique_values(w["extra"]); return chk, n, ctx.expected["unique_values"][w["extra"]]
    raise ValueError(k)


def label(w):
    return w["workload"] if w["workload"] in ("feature", "matrix_sparse") else f"{w['workload']}:{w['extra']}"


def gene_set_index(w):
    return int(w["gene_set"].split(":")[1]) if w["gene_set"] else -1


def run_impl(ctx, impl, opts):
    storage = STORAGE[impl]
    reps = ctx.cfg["reps"]; warm = reps["warm"]; timeout_ms = reps["timeout_s"] * 1000
    rules = ctx.cfg["skip"]["rules"]
    tracemalloc_on = not opts.no_tracemalloc and ctx.manifest["cells"] <= 100_000
    single = bool(opts.single)
    obj, wall, cpu, alloc = timed(lambda: IMPLS[impl](ctx), False)
    if not single:
        ctx.row(implementation=impl, storage=storage, workload="open", phase="open", wall_ms=wall, cpu_ms=cpu, open_ms=wall,
                note=type(obj).__name__)
    all_ids = obj.cells
    works = [parse_single(opts.single)] if single else workloads(ctx)
    seen = set()
    for w in works:
        lab = label(w)
        for r in rules:
            parts = r.split()
            if all((p.startswith("impl=") and p[5:] == impl) or (p.startswith("scale=") and p[6:] == ctx.scale) or
                   (p.startswith("n_features>=") and w["n_features"] >= int(p[12:])) or (p.startswith("workload=") and p[9:] == w["workload"])
                   for p in parts):
                ctx.row(implementation=impl, storage=storage, workload=lab, n_features=w["n_features"], gene_set=gene_set_index(w),
                        cell_subset=w["cell_subset"], layer=w["layer"], status="skipped", note="rule: " + r)
                break
        else:
            family = f"{w['workload']}:{w['layer']}"
            nrep = 1 if single else warm
            rep_start = opts.rep if single else 0
            budget_exceeded = False
            for rep in range(rep_start, rep_start + nrep + (0 if single else 1)):
                phase = ("warmup" if family in seen else "first_call") if rep == 0 else "fetch"
                if not single and rep == 0:
                    seen.add(family)
                key = row_key(impl, lab, w["n_features"], gene_set_index(w), w["cell_subset"], w["layer"], ctx.cache_state, phase, rep)
                if opts.resume and key in ctx.existing:
                    continue
                if budget_exceeded:
                    ctx.row(implementation=impl, storage=storage, workload=lab, n_features=w["n_features"], gene_set=gene_set_index(w),
                            cell_subset=w["cell_subset"], layer=w["layer"], phase=phase, rep=rep, status="skipped_extrapolated",
                            note="first call exceeded timeout budget")
                    continue
                try:
                    res, wall, cpu, alloc = timed(lambda: run_work(ctx, obj, w, all_ids), tracemalloc_on)
                    if res is None:
                        break
                    chk, n, exp = res
                    status, note = "ok", w["extra"]
                except Exception as e:  # noqa: BLE001
                    chk, n, exp, wall, cpu, alloc = "", 0, None, "", "", ""
                    status, note = "error", str(e)[:300]
                ok = checksum_ok(chk, exp) if status == "ok" else ""
                ctx.row(implementation=impl, storage=storage, workload=lab, n_features=w["n_features"], gene_set=gene_set_index(w),
                        cell_subset=w["cell_subset"], layer=w["layer"], n_cells_out=n, phase=phase, rep=rep, wall_ms=wall, cpu_ms=cpu,
                        alloc_mb=alloc, result_checksum=repr(chk) if isinstance(chk, float) else str(chk), checksum_ok=ok, status=status, note=note)
                if status == "ok" and ok == "false":
                    print(f"checksum mismatch: {impl} {lab} {w['gene_set']} {w['cell_subset']} got {chk} expected {exp}", file=sys.stderr)
                if status == "ok" and isinstance(wall, float) and wall > timeout_ms:
                    budget_exceeded = True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scale", required=True)
    ap.add_argument("--config", default=os.path.join(BENCH, "config", "bench.toml"))
    ap.add_argument("--cache", default="warm")
    ap.add_argument("--only", default="")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--run-id", default="")
    ap.add_argument("--single", default="")
    ap.add_argument("--rep", type=int, default=1)
    ap.add_argument("--no-tracemalloc", action="store_true")
    ap.add_argument("--start-time", type=float, default=None)
    opts = ap.parse_args()
    ctx = Ctx(opts)
    impls = opts.only.split(",") if opts.only else list(ctx.cfg["implementations"]["python"])
    t_loaded = time.time()
    for i, impl in enumerate(impls):
        ctx.open_output(f"python_{impl}_{ctx.scale}_{ctx.cache_state}")
        if i == 0 and not opts.single:
            t0 = opts.start_time or T_PROCESS_START
            ctx.row(implementation=impl, storage=STORAGE[impl], workload="startup", phase="startup", rep=0,
                    wall_ms=(t_loaded - t0) * 1000, note="process start -> harness ready")
        print(f"[python {time.strftime('%H:%M:%S')}] {impl} scale={ctx.scale} cache={ctx.cache_state}", flush=True)
        try:
            run_impl(ctx, impl, opts)
        except Exception as e:  # noqa: BLE001
            ctx.row(implementation=impl, storage=STORAGE[impl], workload="open", phase="open", status="error", note=str(e)[:300])
            print(f"  {impl} failed: {e}", file=sys.stderr)
        ctx.close()
    print("results in", ctx.out_dir)


if __name__ == "__main__":
    main()
