# scuba (Python)

`import scuba` gives Python one accessor API over single-cell containers, backed by
[SCUBA.jl](https://github.com/amc-heme/SCUBA.jl) through juliacall.

```python
import anndata as ad, scuba
adata = ad.read_h5ad("AML_h5ad.h5ad")
scuba.fetch_data(adata, ["X_GAPDH", "protein_CD117-AB", "X_umap_1", "Batch"])   # pandas.DataFrame, index = cells
scuba.fetch_data("AML_h5ad.h5ad", ["GAPDH", "nCount_RNA"])                      # native Julia reader, no anndata needed
```

The first call starts Julia and installs SCUBA.jl into an isolated environment (a few minutes
once, seconds afterwards). `import scuba` itself is instant.
