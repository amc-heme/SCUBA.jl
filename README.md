# SCUBA.jl

**S**ingle **C**ell **U**nified **B**ack-end **A**PI, in Julia.

SCUBA.jl gives analysis and visualisation code one accessor API over single-cell containers,
regardless of how they are stored:

```julia
using SCUBA
src = open_h5ad("AML_h5ad.h5ad")
fetch_data(src, ["X_GAPDH", "protein_CD117-AB", "X_umap_1", "Batch"])
```

returns a `DataFrame` with one row per cell and one column per requested variable, whether the
variable is a feature in any modality, a metadata column or a reduction coordinate. The same
call works unchanged on every backend:

| Source | Constructor | Needs |
|---|---|---|
| in-memory container | `MemoryAnnData(...)` | — |
| `.h5ad` (csr / csc / dense, backed) | `open_h5ad(path)` | — |
| AnnData zarr store | `open_zarr(path)` | `using Zarr` |
| BPCells matrix directory | `open_bpcells(dir)` | — |
| live Seurat / SingleCellExperiment object | `rsource(obj)` | `using RCall` |
| live `anndata.AnnData` | `pysource(adata)` | `using PythonCall` |

SCUBA.jl reimplements the R package [SCUBA](https://github.com/amc-heme/SCUBA) (Showers et al.
2024, *F1000Research*) with a single, backend-independent implementation of the variable
resolution rules. It is the engine behind SCUBA for R (v2, through JuliaCall) and the `scuba`
Python package in `python/` (through juliacall).

## Development

```
julia --project=. -e 'using Pkg; Pkg.test()'                                   # core tests
julia --project=. -e 'using Pkg; Pkg.test(test_args=["bridges"])'               # + Zarr / R / Python bridges
```

Bridge tests need R with SeuratObject, SingleCellExperiment and HDF5Array (BPCells optional) and
a Python with anndata, pointed at by `JULIA_PYTHONCALL_EXE` (with `JULIA_CONDAPKG_BACKEND=Null`).
Fixtures are the 250-cell AML CITE-seq dataset in several forms (`fixtures/README.md`); tests look
for them in `fixtures/data`, `$SCUBA_FIXTURES_DIR`, or the lazy `scuba_fixtures` artifact.
Expected outputs recorded from SCUBA R 1.2.2 live in `fixtures/build/goldens`
(`fixtures/make_goldens.R`); `docs/src/compat.md` lists every intentional divergence.
