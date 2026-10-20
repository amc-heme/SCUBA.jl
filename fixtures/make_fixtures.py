#!/usr/bin/env python
"""Generate additional AnnData/MuData forms of the AML fixture from AML_h5ad.h5ad.

Outputs (in the given directory, default fixtures/data):
  AML_h5ad_csr.h5ad    X and layers as CSR
  AML_h5ad_dense.h5ad  X and layers dense
  AML.zarr             zarr store (anndata layout, zarr format 3)
  AML_v2.zarr          zarr store (anndata layout, zarr format 2)
  AML.h5mu             MuData with modalities 'rna' and 'prot'
  AML_h5ad_edge.h5ad   extra obs columns: 'unused_level' (categorical with an unused level),
                       'with_na' (float with NaN), 'nullable_int' (Int64 with <NA>), 'flag' (bool)

Usage: python fixtures/make_fixtures.py [datadir]
Requires: anndata>=0.11, mudata, zarr, scipy, numpy, pandas (conda env `scuba-fixtures`).
"""
import sys, os
import numpy as np, pandas as pd, scipy.sparse as sp
import anndata as ad

datadir = sys.argv[1] if len(sys.argv) > 1 else os.path.join("fixtures", "data")
src = os.path.join(datadir, "AML_h5ad.h5ad")
a = ad.read_h5ad(src)
print(a)

csr = a.copy()
csr.X = sp.csr_matrix(csr.X)
for k in list(csr.layers.keys()):
    csr.layers[k] = sp.csr_matrix(csr.layers[k])
csr.write_h5ad(os.path.join(datadir, "AML_h5ad_csr.h5ad"))

dense = a.copy()
todense = lambda M: np.asarray(M.todense()) if sp.issparse(M) else np.asarray(M)
dense.X = todense(dense.X)
for k in list(dense.layers.keys()):
    dense.layers[k] = todense(dense.layers[k])
dense.write_h5ad(os.path.join(datadir, "AML_h5ad_dense.h5ad"))

try:
    ad.settings.zarr_write_format = 2
    a.write_zarr(os.path.join(datadir, "AML_v2.zarr"))
    ad.settings.zarr_write_format = 3
except Exception as e:
    print("zarr v2 skipped:", e)
a.write_zarr(os.path.join(datadir, "AML.zarr"))

edge = a.copy()
cats = pd.Categorical(edge.obs["Batch"].astype(str), categories=["BM_200AB", "PBMC_200AB", "UNUSED"])
edge.obs["unused_level"] = cats
x = edge.obs["nCount_RNA"].astype(float).to_numpy().copy(); x[::7] = np.nan
edge.obs["with_na"] = x
ni = pd.array(np.arange(edge.n_obs), dtype="Int64"); ni[::5] = pd.NA
edge.obs["nullable_int"] = ni
edge.obs["flag"] = (np.arange(edge.n_obs) % 2 == 0)
edge.write_h5ad(os.path.join(datadir, "AML_h5ad_edge.h5ad"))

try:
    import mudata as md
    rna = ad.AnnData(X=a.X.copy(), obs=a.obs.copy(), var=a.var.copy(), layers={k: v.copy() for k, v in a.layers.items() if k is not None},
                     obsm={"X_pca": a.obsm["X_pca"].copy(), "X_umap": a.obsm["X_umap"].copy()})
    prot_df = a.obsm["protein"]
    prot = ad.AnnData(X=sp.csc_matrix(prot_df.to_numpy(dtype=float)), obs=a.obs[[]].copy(),
                      var=pd.DataFrame(index=prot_df.columns))
    m = md.MuData({"rna": rna, "prot": prot})
    m.obs = a.obs.copy()
    m.write(os.path.join(datadir, "AML.h5mu"))
except Exception as e:  # mudata optional
    print("h5mu skipped:", e)
print("done")
