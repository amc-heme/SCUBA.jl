"""AML_partial.h5mu: a MuData fixture whose modalities hold *different* cell subsets, with a
per-modality obs column that is ambiguous across modalities ("qc"), one unique to a modality
("prot_only"), a modality-level reduction in `prot`, and a global reduction shadowing the bare
alias of `rna`'s X_umap. Built from AML_h5ad.h5ad; run with the `scuba-fixtures` conda env:
    conda run -n scuba-fixtures python fixtures/make_mudata_partial.py
"""
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import mudata as md

datadir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
a = ad.read_h5ad(os.path.join(datadir, "AML_h5ad.h5ad"))
rng = np.random.default_rng(11)

rna = ad.AnnData(X=a.X.copy(), obs=a.obs.copy(), var=a.var.copy(),
                 layers={k: v.copy() for k, v in a.layers.items() if k is not None},
                 obsm={"X_pca": a.obsm["X_pca"].copy(), "X_umap": a.obsm["X_umap"].copy()})
rna = rna[:-3].copy()                                   # 247 cells: the last three are rna-absent
rna.obs["qc"] = rng.normal(size=rna.n_obs)

prot_df = a.obsm["protein"]
prot = ad.AnnData(X=sp.csc_matrix(prot_df.to_numpy(dtype=float)), obs=a.obs[[]].copy(),
                  var=pd.DataFrame(index=prot_df.columns))
keep = np.arange(a.n_obs) % 5 != 0                      # 200 cells: every fifth cell is prot-absent
prot = prot[keep].copy()
prot.obs["qc"] = rng.normal(size=prot.n_obs)            # ambiguous with rna's qc -> keyed only
prot.obs["prot_only"] = pd.Categorical(np.where(np.arange(prot.n_obs) % 2 == 0, "even", "odd"))
prot.obsm["X_pca"] = rng.normal(size=(prot.n_obs, 5))   # modality-level reduction on a partial modality

m = md.MuData({"rna": rna, "prot": prot})
m.obs = a.obs.copy()                                    # global obs over the union (all 250 cells)
m.obsm["X_umap"] = a.obsm["X_umap"].copy()              # global reduction: shadows the bare alias of rna_X_umap
m.write(os.path.join(datadir, "AML_partial.h5mu"))
print("wrote AML_partial.h5mu:", m.shape, {k: v.shape for k, v in m.mod.items()},
      "obs columns:", list(m.obs.columns)[-3:], "obsmap zeros:", {k: int((v == 0).sum()) for k, v in m.obsmap.items()})
