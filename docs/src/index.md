# SCUBA.jl

**S**ingle **C**ell **U**nified **B**ack-end **A**PI for Julia.

SCUBA.jl gives analysis and visualisation code one accessor API over single-cell containers,
regardless of how they are stored. The same call

```julia
using SCUBA
src = open_h5ad("AML_h5ad.h5ad")
fetch_data(src, ["X_GAPDH", "protein_CD117-AB", "X_umap_1", "Batch"])
```

works unchanged on an in-memory `MemoryAnnData`, a backed `.h5ad` file, a zarr store, a Seurat
object with BPCells layers, or a live Seurat / SingleCellExperiment / `anndata.AnnData` object
reached through RCall or PythonCall. It returns a `DataFrame` with one row per cell (the first
column holds the cell IDs) and one column per requested variable, whether that variable is a
feature in any modality, a metadata column or a reduction coordinate.

SCUBA.jl reimplements the R package [SCUBA](https://github.com/amc-heme/SCUBA) (Showers et al.
2024, *F1000Research*) with a single, backend-independent implementation of the variable
resolution rules. It is the engine behind the SCUBA R package (v2) and the `scuba` Python package.

## Installation

```julia
pkg> add SCUBA
```

Optional backends activate when the corresponding package is loaded: `using Zarr` (zarr
stores), `using RCall` (live Seurat / SingleCellExperiment objects), `using PythonCall` (live
`anndata.AnnData` objects), `using Muon` (Muon.jl `AnnData` / `MuData`).

## Contents

```@contents
Pages = ["fetching.md", "backends.md", "compat.md", "api.md"]
Depth = 2
```
