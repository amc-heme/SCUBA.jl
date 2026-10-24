# Backends

| Source | Constructor | Needs | Notes |
|---|---|---|---|
| In-memory container | `MemoryAnnData(...)` | — | layers are cells × features (`SparseMatrixCSC` preferred) or a `LazyMatrix` |
| `.h5ad` file | `open_h5ad(path)` | — | lazy; csr/csc/dense `X` and layers; `obsm` data frames become modalities, arrays become reductions; `raw` becomes modality `raw` |
| zarr store (format 2 or 3) | `open_zarr(path)` | `using Zarr` | same decoder as `.h5ad` |
| BPCells matrix directory | `open_bpcells(dir)` | — | native reader (packed and unpacked); use as a `LazyMatrix` layer or via `SCUBA.bpcells_modality` |
| Live Seurat / SCE object | `rsource(obj)` | `using RCall` | slices in R and scatters them from R memory (no Julia copy); BPCells layers read natively; HDF5Array-backed SCE supported |
| Live `anndata.AnnData` | `pysource(adata)` | `using PythonCall` | in-memory or backed; `set_metadata!` writes through to `adata.obs` |
| MuData `.h5mu` | `open_h5mu(path)` / `open_source(path)` | none (zarr-backed MuData via `using Zarr`) | lazy; modalities `rna_<feature>`, reductions `rna_X_umap` (+ bare alias when unique, global `obsm` arrays by name), metadata = global obs + `<mod>_<col>`; cells absent from a modality read as `NaN` / `missing`; `axis = 0` only |
| Any source + own metadata | `MetadataOverlay(src; obs)` | none | mutable metadata table over a read-only source (`set_metadata!`) |

## Results handed to R (JuliaCall glue)

The R package's Julia engine calls the `ScubaR` glue module inside the RCall extension. Its
`fetch_data` / `fetch_feature` entry points allocate the expression result columns as R double
vectors and fill them in place through the `expression!` sink protocol (`fetch_data(...; alloc)`,
an internal keyword), so a wide result exists once, in R, instead of once in Julia and once in R;
metadata and reduction columns are small and still copied. Cell names never cross the bridge
(the R side already holds them), and `fetch_matrix` returns the sparse pieces of a features × cells
`dgCMatrix`. For 1000 genes × 1M cells this brings the R process peak from ~18.6 GB to ~11–12 GB.

## Aggregation (pseudobulk)

`aggregate(src, codes, ngroups; layer = "counts")` sums (or averages) every feature over groups of
cells. Nothing the size of cells × features is ever built: the request is written into a
`GroupSumSink` (an `ngroups × features` accumulator) through the same `expression!` method every
backend already implements. The cost is therefore the backend's own streaming read:

| Source | How the layer is walked |
|---|---|
| `.h5ad` / zarr csr | runs of consecutive cells, each read bounded by `SCUBA.MAX_READ_ELEMENTS` (all features of a cell arrive together) |
| `.h5ad` / zarr csc | runs of consecutive features, bounded the same way |
| `.h5ad` / zarr dense | feature sub-runs sized so a read stays under `MAX_READ_ELEMENTS` |
| BPCells directory | compressed chunks decoded once, in storage order |
| MuData | the child modality, with cells absent from the modality skipped and not counted (`cell_present`) |
| `rsource` / `pysource` | slices still go through R / Python, so the request is split into cell chunks bounded by `SCUBA.CHUNK_BYTES` (`aggregate_chunk_cells`); native BPCells layers of a Seurat object stream instead |

The R package's `pseudobulk()` calls the `ScubaR.pseudobulk` glue entry point, which returns the
features × groups accumulator as a plain dense R matrix: no sparse index vectors cross the bridge
and the `Int32` nonzero limit of `fetch_matrix` does not apply.

## Writing a backend

Subtype `AbstractSingleCellSource` and implement the required methods listed in its docstring.
The core never asks for whole matrices: `expression(src, modality, layer, feature_idx, cell_idx)`
receives integer indices for exactly the block it needs. For speed and memory, implement
`expression!(out, src, modality, layer, feature_idx, cell_idx)` instead and write every value
with `SCUBA._put!(out, i, j, v)`: `out` is a dense matrix, a vector of per-feature columns or a
`SparseSink`, so one method yields typed dense, column and sparse results (`expression`,
`expression_columns` and `expression_sparse` then come for free; override `expression_sparse`
only with a cheaper native sparse path). Other optional methods (`supports_sparse`,
`metadata_table`, `set_metadata!`, `canonical_modality`, ...) have defaults.
Run `test_source_conformance` from `test/conformance.jl` against a new backend.
