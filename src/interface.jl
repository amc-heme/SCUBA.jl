# The backend protocol. Every container type that SCUBA can read implements the
# required methods below; the accessor API (api.jl) is written once on top of them.

"""
    AbstractSingleCellSource

Supertype of every SCUBA backend. Required methods:

- `cell_ids(src)`, `modality_keys(src)`, `default_modality(src)`, `default_layer(src)`
- `feature_names(src, modality)`, `layer_names(src, modality)`
- `expression(src, modality, layer, feature_idx, cell_idx)` → `Matrix{Float64}` (cells × features)
- `metadata_names(src)`, `metadata_column(src, name, cell_idx)`
- `reduction_keys(src)`, `reduction_ndims(src, key)`, `reduction(src, key, dims, cell_idx)` → `Matrix{Float64}` (cells × dims)

Optional methods have defaults: `supports_sparse`, `expression_sparse`, `reduction_dimnames`,
`layer_exists`, `metadata_table`, `set_metadata!`, `invalidate!`, `cache`, `source_name`.
"""
abstract type AbstractSingleCellSource end

# ---------------------------------------------------------------------------
# Per-source cache of name → index dictionaries

mutable struct SourceCache
    cell_index::Union{Nothing,Dict{String,Int}}
    feature_index::Dict{String,Dict{String,Int}}
    extra::Dict{Symbol,Any}
    lock::ReentrantLock
end
SourceCache() = SourceCache(nothing, Dict{String,Dict{String,Int}}(), Dict{Symbol,Any}(), ReentrantLock())

function Base.empty!(c::SourceCache)
    lock(c.lock) do
        c.cell_index = nothing
        empty!(c.feature_index)
        empty!(c.extra)
    end
    return c
end

_notimpl(method, src) = throw(NotImplementedBySource(String(method), string(nameof(typeof(src)))))

# ---------------------------------------------------------------------------
# Required interface

"""`cell_ids(src)` — cell identifiers in object order."""
cell_ids(src::AbstractSingleCellSource) = _notimpl(:cell_ids, src)
"""`modality_keys(src)` — modality keys, main modality first, without a trailing underscore."""
modality_keys(src::AbstractSingleCellSource) = _notimpl(:modality_keys, src)
default_modality(src::AbstractSingleCellSource) = first(modality_keys(src))
"""`default_layer(src)` — the normalised-data layer name (Seurat `data`, SCE `logcounts`, AnnData `X`)."""
default_layer(src::AbstractSingleCellSource) = _notimpl(:default_layer, src)
feature_names(src::AbstractSingleCellSource, modality::AbstractString) = _notimpl(:feature_names, src)
layer_names(src::AbstractSingleCellSource, modality::AbstractString) = _notimpl(:layer_names, src)
"""
    expression(src, modality, layer, feature_idx, cell_idx) -> Matrix{Float64}

Dense block of expression values, cells as rows (in `cell_idx` order) and features as
columns (in `feature_idx` order).
"""
expression(src::AbstractSingleCellSource, modality::AbstractString, layer::AbstractString,
           feature_idx::AbstractVector{Int}, cell_idx::AbstractVector{Int}) = _notimpl(:expression, src)
"""
    expression_columns(src, modality, layer, feature_idx, cell_idx, T) -> Vector{Vector{T}}

Per-feature result columns (each of length `length(cell_idx)`, in `cell_idx` order). The default
materialises `expression` and copies; backends override it to fill the columns directly so a
large request never holds a second dense copy.
"""
function expression_columns(src::AbstractSingleCellSource, modality::AbstractString, layer::AbstractString,
                            feature_idx::AbstractVector{Int}, cell_idx::AbstractVector{Int}, ::Type{T}) where {T<:Real}
    M = expression(src, modality, layer, feature_idx, cell_idx)
    return Vector{T}[Vector{T}(view(M, :, j)) for j in 1:size(M, 2)]
end
metadata_names(src::AbstractSingleCellSource) = _notimpl(:metadata_names, src)
"""`metadata_column(src, name, cell_idx)` — one metadata column, subset to `cell_idx`."""
metadata_column(src::AbstractSingleCellSource, name::AbstractString, cell_idx::AbstractVector{Int}) =
    _notimpl(:metadata_column, src)
reduction_keys(src::AbstractSingleCellSource) = _notimpl(:reduction_keys, src)
reduction_ndims(src::AbstractSingleCellSource, key::AbstractString) = _notimpl(:reduction_ndims, src)
"""`reduction(src, key, dims, cell_idx)` — embedding block, cells × dims."""
reduction(src::AbstractSingleCellSource, key::AbstractString, dims::AbstractVector{Int},
          cell_idx::AbstractVector{Int}) = _notimpl(:reduction, src)

# ---------------------------------------------------------------------------
# Optional interface

cache(src::AbstractSingleCellSource) = src.cache::SourceCache
source_name(src::AbstractSingleCellSource) = string(nameof(typeof(src)))
supports_sparse(src::AbstractSingleCellSource, modality::AbstractString, layer::AbstractString) = false
"""`expression_sparse(src, modality, layer, feature_idx, cell_idx)` — as `expression`, but a `SparseMatrixCSC`."""
expression_sparse(src::AbstractSingleCellSource, modality::AbstractString, layer::AbstractString,
                  feature_idx::AbstractVector{Int}, cell_idx::AbstractVector{Int}) =
    sparse(expression(src, modality, layer, feature_idx, cell_idx))
"""`reduction_search_names(src, key)` — names under which [`default_reduction`](@ref) recognises a reduction (default: the key only; the Seurat bridge adds the reduction name)."""
reduction_search_names(src::AbstractSingleCellSource, key::AbstractString) = String[String(key)]

"""`reduction_dimnames_impl(src, key, dims)` — backend hook naming the dimensions of a reduction (default `key_dim`)."""
reduction_dimnames_impl(src::AbstractSingleCellSource, key::AbstractString, dims::AbstractVector{Int}) =
    [string(key, "_", d) for d in dims]
layer_exists(src::AbstractSingleCellSource, modality::AbstractString, layer::AbstractString) =
    layer in layer_names(src, modality)
"""`set_metadata!(src, table)` — replace the metadata table. Bridges write through to the host object."""
set_metadata!(src::AbstractSingleCellSource, table) = throw(ReadOnlySource(source_name(src)))
"""`invalidate!(src)` — drop cached name indices after the underlying object was mutated."""
invalidate!(src::AbstractSingleCellSource) = (empty!(cache(src)); src)

function metadata_table(src::AbstractSingleCellSource, cell_idx::AbstractVector{Int})
    df = DataFrame()
    for name in metadata_names(src)
        df[!, name] = metadata_column(src, name, cell_idx)
    end
    return df
end

# ---------------------------------------------------------------------------
# Derived helpers (cached)

ncells(src::AbstractSingleCellSource) = length(cell_ids(src))

function cell_index(src::AbstractSingleCellSource)
    c = cache(src)
    c.cell_index === nothing || return c.cell_index
    lock(c.lock) do
        if c.cell_index === nothing
            ids = cell_ids(src)
            d = Dict{String,Int}()
            sizehint!(d, length(ids))
            for (i, id) in enumerate(ids)
                d[String(id)] = i   # first occurrence wins for duplicated IDs
            end
            c.cell_index = d
        end
    end
    return c.cell_index
end

function feature_index(src::AbstractSingleCellSource, modality::AbstractString)
    c = cache(src)
    d = get(c.feature_index, modality, nothing)
    d === nothing || return d
    lock(c.lock) do
        get!(c.feature_index, String(modality)) do
            names = feature_names(src, modality)
            d = Dict{String,Int}()
            sizehint!(d, length(names))
            for (i, n) in enumerate(names)
                get!(d, String(n), i)
            end
            d
        end
    end
end

"""
    canonical_modality(src, name) -> String

Map an alternative modality name (e.g. a Seurat assay name `"RNA"` for key `"rna"`) to the
modality key. The default is the identity; backends with aliases override it.
"""
canonical_modality(src::AbstractSingleCellSource, name::AbstractString) = String(name)
"""`canonical_reduction(src, name)` — as [`canonical_modality`](@ref) for reductions (e.g. Seurat `"umap"` → key `"UMAP"`)."""
canonical_reduction(src::AbstractSingleCellSource, name::AbstractString) = String(name)

function check_modality(src::AbstractSingleCellSource, modality::AbstractString)
    mods = modality_keys(src)
    key = canonical_modality(src, modality)
    key in mods || throw(ModalityNotFound(String(modality), collect(mods)))
    return key
end

function check_reduction(src::AbstractSingleCellSource, red::AbstractString)
    reds = reduction_keys(src)
    key = canonical_reduction(src, red)
    key in reds || throw(ReductionNotFound(String(red), collect(reds)))
    return key
end

"""
    resolve_layer(src, modality, layer) -> String

Return the layer to read for `modality`. `nothing` means the source default. If the
requested layer does not exist but the modality has exactly one layer, that layer is used
with a `:layer_substituted` warning (modalities such as AnnData `obsm` data frames have no
layers of their own).
"""
function resolve_layer(src::AbstractSingleCellSource, modality::AbstractString, layer)
    layers = layer_names(src, modality)
    requested = layer === nothing ? default_layer(src) : String(layer)
    requested in layers && return requested
    if length(layers) == 1
        if layer !== nothing
            scuba_warn(:layer_substituted,
                       "Layer '$requested' does not exist in modality '$modality'; using its only layer '$(layers[1])'.")
        end
        return String(layers[1])
    end
    throw(LayerNotFound(requested, String(modality), collect(layers)))
end

function Base.show(io::IO, ::MIME"text/plain", src::AbstractSingleCellSource)
    print(io, source_name(src), " with ", ncells(src), " cells")
    mods = modality_keys(src)
    println(io)
    for m in mods
        print(io, "  modality '", m, "': ", length(feature_names(src, m)), " features, layers ",
              join(layer_names(src, m), ", "))
        m == default_modality(src) && print(io, " (default)")
        println(io)
    end
    reds = reduction_keys(src)
    isempty(reds) || println(io, "  reductions: ", join(["$r ($(reduction_ndims(src, r)))" for r in reds], ", "))
    print(io, "  metadata: ", length(metadata_names(src)), " columns")
end

# ---------------------------------------------------------------------------
# Constructors provided by package extensions (documented here; methods added when the
# corresponding package is loaded)

"""
    rsource(obj)

Wrap a live R Seurat or SingleCellExperiment object (an `RCall.RObject`, or a `Symbol` naming
a variable in the R session) as a SCUBA source. Requires `using RCall`.
"""
function rsource end

"""
    pysource(obj)

Wrap a live Python `anndata.AnnData` object (a `PythonCall.Py`) as a SCUBA source. Requires
`using PythonCall`.
"""
function pysource end

"""
    open_zarr(path)

Open an AnnData zarr store as a SCUBA source. Requires `using Zarr`.
"""
function open_zarr end

"""
    to_pandas(df)

Convert a SCUBA result `DataFrame` to a pandas DataFrame (cell IDs as index). Requires `using PythonCall`.
"""
function to_pandas end
