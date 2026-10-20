#!/usr/bin/env python
"""Generate one synthetic single-cell dataset for the SCUBA benchmark suite.

Writes, under <data_dir>/<scale>/:
  csc.h5ad, csr.h5ad      X = log1p(CP10k) float32, layers/counts float32, obs (20 columns), obsm X_umap / X_pca
  obs.parquet, obsm_X_umap.parquet, obsm_X_pca.parquet   (for the R object builder)
  selection.json          seeded gene sets and cell subsets shared by every harness
  expected.json           reference checksums for every workload
  MANIFEST.json           parameters, file sizes, timings; `complete: true` marks a finished dataset

Usage: python -I generate/make_dataset.py --scale 10k [--config config/bench.toml] [--data-dir DIR] [--force]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
import tomllib

import numpy as np
import pandas as pd
import scipy.sparse as sp

GENERATOR_VERSION = "0.1.0"
CELL_TYPES = ["T_cell", "B_cell", "NK_cell", "CD14_Monocyte", "CD16_Monocyte", "cDC", "pDC", "Plasma_cell",
              "HSC", "Erythroid", "Megakaryocyte", "Neutrophil", "Basophil", "Fibroblast", "Endothelial"]
CELL_TYPE_P = np.array([0.25, 0.12, 0.08, 0.15, 0.05, 0.03, 0.02, 0.03, 0.04, 0.07, 0.03, 0.06, 0.02, 0.03, 0.02])


def log(msg):
    print(f"[make_dataset {time.strftime('%H:%M:%S')}] {msg}", flush=True)


def sha1_lines(values) -> str:
    return hashlib.sha1("\n".join(map(str, values)).encode()).hexdigest()


def git_sha(path):
    try:
        return subprocess.run(["git", "-C", path, "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return None


def generate_blocks(rng, n_cells, n_genes, gene_mean, nnz_target, block):
    """Yield (indices int32, counts float32, k) per block; indices sorted within each cell."""
    p = gene_mean / gene_mean.sum()
    gene_logp = np.log(p).astype(np.float32)
    start = 0
    while start < n_cells:
        b = min(block, n_cells - start)
        k = rng.poisson(nnz_target, b).clip(50, n_genes).astype(np.int64)
        kmax = int(k.max())
        size = np.exp(rng.normal(0.0, SIZE_SIGMA, b)).astype(np.float32)
        g = rng.gumbel(size=(b, n_genes)).astype(np.float32)
        g += gene_logp[None, :]
        cand = np.argpartition(-g, kmax - 1, axis=1)[:, :kmax]              # unordered top-kmax
        cand_scores = np.take_along_axis(g, cand, axis=1)
        del g
        order = np.argsort(-cand_scores, axis=1)
        cand = np.take_along_axis(cand, order, axis=1)                       # descending score
        del cand_scores, order
        valid = np.arange(kmax)[None, :] < k[:, None]
        cand = np.where(valid, cand, n_genes)                                # sentinel sorts last
        cand.sort(axis=1)
        idx = cand[valid.sum(axis=1)[:, None] > np.arange(kmax)[None, :]]    # row-major flatten of valid entries
        # idx is already row-major because `valid` counts are prefixes after the sort
        lam = gene_mean[idx] * np.repeat(size, k)
        counts = (1 + rng.poisson(lam)).astype(np.float32)
        yield idx.astype(np.int32), counts, k
        start += b
        log(f"  cells {start:,}/{n_cells:,}")


SIZE_SIGMA = 0.3


def build_obs(rng, n_cells, counts_rowsum, nfeat):
    ct = rng.choice(len(CELL_TYPES), n_cells, p=CELL_TYPE_P / CELL_TYPE_P.sum())
    df = pd.DataFrame(index=pd.Index([f"cell_{i:010d}" for i in range(1, n_cells + 1)], name="cell_id"))
    df["cell_type"] = pd.Categorical.from_codes(ct, CELL_TYPES)
    df["batch"] = pd.Categorical.from_codes(rng.integers(0, 8, n_cells), [f"b{i}" for i in range(1, 9)])
    df["sample"] = pd.Categorical.from_codes(rng.integers(0, 48, n_cells), [f"s{i:02d}" for i in range(1, 49)])
    df["donor"] = pd.Categorical.from_codes(rng.integers(0, 12, n_cells), [f"d{i:02d}" for i in range(1, 13)])
    df["condition"] = pd.Categorical.from_codes(rng.integers(0, 3, n_cells), ["ctrl", "treatA", "treatB"])
    df["cluster"] = pd.Categorical.from_codes(rng.integers(0, 30, n_cells), [f"c{i}" for i in range(30)])
    df["phase"] = pd.Categorical.from_codes(rng.integers(0, 3, n_cells), ["G1", "S", "G2M"])
    df["nCount_RNA"] = counts_rowsum.astype(np.float64)
    df["nFeature_RNA"] = nfeat.astype(np.int64)
    df["percent_mt"] = (rng.beta(2, 30, n_cells) * 100).astype(np.float64)
    df["doublet_score"] = rng.beta(1, 20, n_cells).astype(np.float64)
    df["score_a"] = rng.normal(size=n_cells)
    df["score_b"] = rng.normal(size=n_cells)
    df["score_c"] = rng.uniform(size=n_cells)
    df["n_genes_by_counts"] = nfeat.astype(np.int64)
    df["total_counts_mt"] = rng.integers(0, 500, n_cells).astype(np.int64)
    df["leiden_int"] = rng.integers(0, 40, n_cells).astype(np.int64)
    df["is_doublet"] = df["doublet_score"].to_numpy() > 0.3
    df["pass_qc"] = df["percent_mt"].to_numpy() < 15
    df["barcode"] = [f"{rng.integers(0, 2**32):08X}-1" for _ in range(min(n_cells, 1))] * 0 + \
        np.char.add(np.char.add(np.array([f"{x:08X}" for x in rng.integers(0, 2**32, n_cells)]), "-"),
                    (ct + 1).astype(str)).tolist()
    return df, ct


def build_obsm(rng, n_cells, ct):
    angles = np.linspace(0, 2 * np.pi, len(CELL_TYPES), endpoint=False)
    centers = np.stack([10 * np.cos(angles), 10 * np.sin(angles)], axis=1)
    umap = (centers[ct] + rng.normal(0, 1.2, (n_cells, 2))).astype(np.float32)
    pca = rng.normal(0, 1, (n_cells, 50)).astype(np.float32)
    pca[:, 0] += centers[ct, 0] / 2
    pca[:, 1] += centers[ct, 1] / 2
    return umap, pca


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scale", required=True)
    ap.add_argument("--config", default=os.path.join(os.path.dirname(__file__), "..", "config", "bench.toml"))
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    with open(args.config, "rb") as fh:
        cfg = tomllib.load(fh)
    sc = cfg["scales"][args.scale]
    gen = cfg["generation"]
    wl = cfg["workloads"]
    data_root = os.path.expanduser(args.data_dir or cfg["suite"]["data_dir"])
    out = os.path.join(data_root, args.scale)
    os.makedirs(out, exist_ok=True)
    manifest_path = os.path.join(out, "MANIFEST.json")
    if os.path.exists(manifest_path) and not args.force:
        with open(manifest_path) as fh:
            if json.load(fh).get("complete"):
                log(f"{args.scale}: dataset complete, skipping (use --force to regenerate)")
                return

    n_cells, n_genes, nnz_target, seed = sc["cells"], sc["genes"], sc["nnz_per_cell"], sc["seed"]
    global SIZE_SIGMA
    SIZE_SIGMA = gen["size_factor_log_sigma"]
    rng = np.random.default_rng(seed)
    t_all = time.time()
    timings = {}

    gene_names = np.array([f"gene{i:05d}" for i in range(1, n_genes + 1)])  # no underscores: Seurat rewrites them
    gene_mean = np.exp(rng.normal(gen["gene_mean_log_mu"], gen["gene_mean_log_sigma"], n_genes)).astype(np.float64)

    log(f"{args.scale}: generating {n_cells:,} cells x {n_genes:,} genes, ~{nnz_target} nnz/cell")
    t0 = time.time()
    idx_parts, cnt_parts, k_parts = [], [], []
    for idx, cnt, k in generate_blocks(rng, n_cells, n_genes, gene_mean, nnz_target, gen["block_cells"]):
        idx_parts.append(idx); cnt_parts.append(cnt); k_parts.append(k)
    k_all = np.concatenate(k_parts)
    indptr = np.zeros(n_cells + 1, dtype=np.int64)
    np.cumsum(k_all, out=indptr[1:])
    indices = np.concatenate(idx_parts); del idx_parts
    counts = np.concatenate(cnt_parts); del cnt_parts
    nnz = int(indptr[-1])
    timings["generate_s"] = time.time() - t0
    log(f"nnz = {nnz:,} ({timings['generate_s']:.1f}s)")

    # per-cell totals and log-normalised layer sharing the sparsity pattern
    rowsum = np.add.reduceat(counts.astype(np.float64), indptr[:-1])
    scale_per_cell = (1e4 / rowsum).astype(np.float32)
    xdata = np.log1p(counts * np.repeat(scale_per_cell, k_all)).astype(np.float32)

    csr_counts = sp.csr_matrix((counts, indices, indptr), shape=(n_cells, n_genes))
    csr_x = sp.csr_matrix((xdata, indices.copy(), indptr.copy()), shape=(n_cells, n_genes))

    obs, ct = build_obs(rng, n_cells, rowsum, k_all)
    umap, pca = build_obsm(rng, n_cells, ct)
    cell_ids = obs.index.to_numpy()

    import anndata as ad
    # plain string arrays and natural index dtypes keep the files readable by BPCells, HDF5Array and SCUBA.jl
    ad.settings.allow_write_nullable_strings = False
    ad.settings.write_csr_csc_indices_with_min_possible_dtype = False
    var = pd.DataFrame(index=pd.Index(gene_names, name="gene"))
    var["gene_mean"] = gene_mean

    t0 = time.time()
    csr_layers = sc.get("csr_layers", ["X", "counts"])
    layers = {"counts": csr_counts} if "counts" in csr_layers else {}
    adata = ad.AnnData(X=csr_x, obs=obs, var=var, layers=layers, obsm={"X_umap": umap, "X_pca": pca})
    adata.write_h5ad(os.path.join(out, "csr.h5ad"), compression=None)
    del adata, layers
    timings["write_csr_s"] = time.time() - t0
    log(f"wrote csr.h5ad ({timings['write_csr_s']:.1f}s)")

    t0 = time.time()
    csc_x = csr_x.tocsc(); del csr_x
    csc_counts = csr_counts.tocsc(); del csr_counts
    adata = ad.AnnData(X=csc_x, obs=obs, var=var, layers={"counts": csc_counts}, obsm={"X_umap": umap, "X_pca": pca})
    adata.write_h5ad(os.path.join(out, "csc.h5ad"), compression=None)
    del adata
    timings["write_csc_s"] = time.time() - t0
    log(f"wrote csc.h5ad ({timings['write_csc_s']:.1f}s)")

    # tables for the R object builder
    obs.reset_index().to_parquet(os.path.join(out, "obs.parquet"), index=False)
    pd.DataFrame(umap, columns=["UMAP_1", "UMAP_2"]).to_parquet(os.path.join(out, "obsm_X_umap.parquet"), index=False)
    pd.DataFrame(pca, columns=[f"PC_{i}" for i in range(1, 51)]).to_parquet(os.path.join(out, "obsm_X_pca.parquet"), index=False)

    # deterministic selections
    g_rng = np.random.default_rng(gen["gene_selection_seed"])
    gene_sets = {}
    for n in wl["n_features"]:
        for rep in range(wl["gene_set_reps"]):
            gene_sets[f"{n}:{rep}"] = g_rng.choice(gene_names, n, replace=False).tolist()
    c_rng = np.random.default_rng(seed + 1)
    subsets = {
        "all": None,
        "rand1": np.sort(c_rng.choice(n_cells, max(1, n_cells // 100), replace=False)).tolist(),
        "rand10": np.sort(c_rng.choice(n_cells, max(1, n_cells // 10), replace=False)).tolist(),
        "contig1": [0, max(1, n_cells // 100)],
        "contig10": [0, max(1, n_cells // 10)],
    }
    selection = {"gene_sets": gene_sets, "cell_subsets": subsets, "n_cells": n_cells, "n_genes": n_genes,
                 "note": "cell subset indices are 0-based; contig* entries are [start, stop)"}
    with open(os.path.join(out, "selection.json"), "w") as fh:
        json.dump(selection, fh)

    # expected checksums
    t0 = time.time()
    gene_pos = {g: i for i, g in enumerate(gene_names)}

    def cell_index(name):
        s = subsets[name]
        if s is None:
            return None
        if name.startswith("contig"):
            return np.arange(s[0], s[1])
        return np.asarray(s)

    expected = {"feature": {}, "metadata": {}, "reduction": {}, "fetch_cells": {}, "unique_values": {}}
    for layer_name, M in (("data", csc_x), ("counts", csc_counts)):
        for key, genes in gene_sets.items():
            n, rep = key.split(":")
            if layer_name == "counts" and rep != "0":
                continue
            gidx = np.array([gene_pos[g] for g in genes])
            cols = M[:, gidx].tocsr()
            for sub in wl["cell_subsets"]:
                ci = cell_index(sub)
                block = cols if ci is None else cols[ci]
                total = float(block.sum(dtype=np.float64))
                coo = block.tocoo()
                if coo.nnz:
                    order = np.lexsort((coo.col, coo.row))
                    j = order[0]
                    first = [str(cell_ids[ci[coo.row[j]] if ci is not None else coo.row[j]]), str(genes[coo.col[j]]), float(coo.data[j])]
                else:
                    first = None
                expected["feature"][f"{layer_name}|{key}|{sub}"] = {"sum": total, "nnz": int(block.nnz), "first": first,
                                                                    "n_cells": int(block.shape[0])}
    expected["metadata"]["full"] = f"{n_cells}:{float(obs['nCount_RNA'].sum())!r}"   # ncol differs per class (Seurat adds orig.ident)
    expected["metadata"]["nCount_RNA"] = f"{n_cells}:{float(obs['nCount_RNA'].sum())!r}"
    expected["metadata"]["cell_type"] = f"{n_cells}:{sha1_lines(obs['cell_type'].astype(str).tolist())}"
    expected["reduction"]["X_umap|1-2"] = float(umap.astype(np.float64).sum())
    expected["reduction"]["X_pca|1-2"] = float(pca[:, :2].astype(np.float64).sum())
    expected["reduction"]["X_pca|1-50"] = float(pca.astype(np.float64).sum())
    for spec in wl["fetch_cells"]:
        col, level = spec.split(":")
        ids = sorted(cell_ids[(obs[col].astype(str) == level).to_numpy()].tolist())
        expected["fetch_cells"][spec] = f"{len(ids)}:{sha1_lines(ids)}"
    for col in wl["unique_values"]:
        lv = sorted(set(obs[col].astype(str).tolist()))
        expected["unique_values"][col] = f"{len(lv)}:{sha1_lines(lv)}"
    with open(os.path.join(out, "expected.json"), "w") as fh:
        json.dump(expected, fh)
    timings["expected_s"] = time.time() - t0

    files = {f: os.path.getsize(os.path.join(out, f)) for f in sorted(os.listdir(out)) if f != "MANIFEST.json"}
    manifest = {
        "scale": args.scale, "cells": n_cells, "genes": n_genes, "nnz": nnz, "nnz_per_cell": nnz / n_cells,
        "seed": seed, "generator_version": GENERATOR_VERSION, "generation": gen, "csr_layers": csr_layers,
        "files": files, "timings_s": timings, "total_s": time.time() - t_all,
        "python": sys.version.split()[0], "anndata": __import__("importlib.metadata").metadata.version("anndata"),
        "numpy": np.__version__, "scipy": sp.__version__ if hasattr(sp, "__version__") else None,
        "platform": platform.platform(), "scuba_jl_sha": git_sha(os.path.join(os.path.dirname(__file__), "..", "..")),
        "created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "complete": True,
    }
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    log(f"done in {manifest['total_s']:.1f}s; {sum(files.values()) / 2**30:.2f} GiB on disk")


if __name__ == "__main__":
    main()
