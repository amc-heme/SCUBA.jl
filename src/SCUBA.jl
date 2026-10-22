"""
SCUBA.jl — Single Cell Unified Back-end API.

One accessor API (`fetch_data`, `fetch_metadata`, `fetch_reduction`, ...) over single-cell
containers stored as AnnData (in memory, .h5ad, zarr), BPCells directories, or live R and
Python objects through package extensions.
"""
module SCUBA

using DataFrames
using CategoricalArrays
using Tables
using OrderedCollections
using SparseArrays
using Logging
using HDF5
using Mmap
using Artifacts, LazyArtifacts

include("errors.jl")
include("interface.jl")
include("keys.jl")
include("resolve.jl")
include("sparse.jl")
include("table.jl")
include("api.jl")
include("chunks.jl")
include("memory.jl")
include("overlay.jl")
include("encoding/store.jl")
include("encoding/anndata_spec.jl")
include("h5ad.jl")
include("mudata.jl")
include("open.jl")
include("bpcells/bp128.jl")
include("bpcells/matrix.jl")
include("fixtures.jl")
include("precompile.jl")

# core API
export fetch_data, fetch_feature, fetch_matrix, fetch_metadata, fetch_reduction, fetch_cells,
       all_cells, get_all_cells, meta_varnames, all_keys, features_in_modality, features_in_assay,
       reduction_dimnames, default_layer, default_reduction, default_modality, unique_values, n_unique,
       set_metadata!, update_object_metadata!, snapshot, resolve_vars, cell_column, cells, drop_cells, eachchunk
# sources
export AbstractSingleCellSource, MemoryAnnData, Modality, MetadataOverlay, H5ADSource, open_h5ad, MuDataSource, open_h5mu, open_source, BPCellsMatrix, open_bpcells,
       rsource, pysource, open_zarr, to_pandas, to_scipy, to_numpy
# errors & warnings
export SCUBAError, VarsNotFound, MetadataNotFound, LayerNotFound, ModalityNotFound, ReductionNotFound,
       NoDefaultReduction, InvalidDims, CellsNotFound, NonCategoricalVar, ReadOnlySource, ResultTooLarge,
       SCUBAWarning, with_collected_warnings
# fixtures
export fixtures_dir, fixture_path

end # module
