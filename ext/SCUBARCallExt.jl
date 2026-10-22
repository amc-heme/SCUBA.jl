# Live R objects (Seurat, SingleCellExperiment) as SCUBA sources, through RCall.
#
# Nothing is copied at construction except names and shapes. Expression blocks are sliced in R
# and only the requested (features × cells) rectangle crosses into Julia. BPCells layers that are
# plain on-disk matrices are read natively by `BPCellsMatrix`, bypassing R entirely.
module SCUBARCallExt

using SCUBA, RCall, DataFrames, CategoricalArrays, SparseArrays, OrderedCollections
using SCUBA: AbstractSingleCellSource, SourceCache, LazyMatrix, BPCellsMatrix, pull, pull!, pull_sparse, pull_columns,
             gather_features_by_cells, gather_fc_sparse, gather_fc!, _columns, _put!, ModalityNotFound, LayerNotFound,
             ReductionNotFound, ReadOnlySource
import SCUBA: cell_ids, modality_keys, default_modality, default_layer, feature_names, layer_names,
              expression, expression!, expression_columns, expression_sparse, supports_sparse, metadata_names, metadata_column, metadata_table,
              reduction_keys, reduction_ndims, reduction, source_name, set_metadata!, invalidate!,
              canonical_modality, canonical_reduction, reduction_search_names

const R_HELPERS = """
scuba_bp_dir <- function(x) {
  # Return list(dir, transpose) for a BPCells matrix that is a plain MatrixDir, possibly wrapped
  # only in dimension-renaming wrappers; NULL otherwise (transformed/queued matrices need R).
  # BPCells implements t() by flipping the `transpose` slot (set consistently on every level of a
  # wrapper chain). The flag alone does not say whether the files are transposed relative to the
  # logical matrix: directories written by transpose_storage_order() carry storage_order "row"
  # and transpose = TRUE while their shape file is the logical shape. The Julia side resolves the
  # orientation from the layer's dimensions and uses this flag only as a tie-breaker.
  transpose <- methods::.hasSlot(x, "transpose") && isTRUE(x@transpose)
  while (!methods::is(x, "MatrixDir")) {
    cls <- class(x)[1]
    if (!(cls %in% c("RenameDims")) || !methods::.hasSlot(x, "matrix")) return(NULL)
    x <- x@matrix
  }
  if (!methods::.hasSlot(x, "dir")) return(NULL)
  list(dir = x@dir, transpose = transpose)
}
scuba_seurat_layers <- function(assay) {
  if (methods::is(assay, "Assay5")) return(SeuratObject::Layers(assay))
  slots <- c("counts", "data", "scale.data")
  slots[vapply(slots, function(s) {
    m <- tryCatch(SeuratObject::GetAssayData(assay, layer = s), error = function(e) NULL)
    !is.null(m) && length(m) > 0 && nrow(m) > 0
  }, logical(1))]
}
scuba_seurat_layer_data <- function(assay, layer) {
  if (methods::is(assay, "Assay5")) SeuratObject::LayerData(assay, layer = layer)
  else SeuratObject::GetAssayData(assay, layer = layer)
}
scuba_slice <- function(m, f, c) {
  s <- m[f, c, drop = FALSE]
  if (methods::is(s, "DelayedArray") || methods::is(s, "IterableMatrix")) s <- as.matrix(s)
  if (!methods::is(s, "dgCMatrix") && !is.matrix(s)) s <- as.matrix(s)
  s
}
"""

function __init__()
    reval(R_HELPERS)
    return nothing
end

_rstrings(x) = x === nothing || isa(x, Nothing) ? String[] : collect(String, rcopy(Vector{String}, x))
_rint(v::AbstractVector{Int}) = RCall.robject(Int32.(v))

# Convert an R matrix / dgCMatrix slice (features × cells) to the cells × features block.
_is_dgc(s::RObject) = rcopy(Bool, rcall(:inherits, s, "dgCMatrix"))
# dgCMatrix (features × cells) -> SparseMatrixCSC without relying on RCall sparse conversions
function _dgc_to_sparse(s::RObject)
    i = rcopy(Vector{Int32}, R"$s@i") .+ 1
    p = rcopy(Vector{Int32}, R"$s@p") .+ 1
    x = rcopy(Vector{Float64}, R"as.numeric($s@x)")
    d = rcopy(Vector{Int}, R"$s@Dim")
    return SparseMatrixCSC{Float64,Int}(d[1], d[2], Vector{Int}(p), Vector{Int}(i), x)
end
# R returns a features × cells slice as a matrix, or as a bare numeric when it is 1 × 1 / a
# vector; reshape by the known dimensions (R matrices are column-major like Julia's).
_dense_from_r(s::RObject, nf::Int, nc::Int, ::Type{T}) where {T} = reshape(rcopy(Vector{T}, s), nf, nc)
function _block_from_r(s::RObject, nf::Int, nc::Int)
    if _is_dgc(s)
        return gather_features_by_cells(_dgc_to_sparse(s), collect(1:nc), collect(1:nf))
    else
        return Matrix{Float64}(permutedims(_dense_from_r(s, nf, nc, Float64)))
    end
end
# Typed per-feature columns from an R slice (features × cells); the R-sliced values are copied into
# `T` directly so a large request never holds a Float64 intermediate on the Julia side.
function _columns_from_r(s::RObject, nf::Int, nc::Int, ::Type{T}) where {T}
    cols = _columns(T, nc, nf)
    if _is_dgc(s)
        gather_fc!(cols, _dgc_to_sparse(s), collect(1:nc), collect(1:nf))
    else
        M = _dense_from_r(s, nf, nc, T)               # features × cells, converted while copying
        for j in 1:nf
            copyto!(cols[j], view(M, j, :))
        end
    end
    return cols
end
function _sparse_block_from_r(s::RObject, nf::Int, nc::Int)
    _is_dgc(s) && return gather_fc_sparse(_dgc_to_sparse(s), collect(1:nc), collect(1:nf))
    return sparse(Matrix{Float64}(permutedims(_dense_from_r(s, nf, nc, Float64))))
end
# Scatter an R slice (features × cells) into any `_put!` sink straight from R memory: the dgCMatrix
# slots / dense vector are wrapped with `unsafe_vec`, so no Julia copy of the slice is made.
function _scatter_from_r!(out, s::RObject, nf::Int, nc::Int)
    if _is_dgc(s)
        ri = R"$s@i"; rp = R"$s@p"; rx = R"$s@x"
        GC.@preserve ri rp rx begin
            i = RCall.unsafe_vec(ri); p = RCall.unsafe_vec(rp); x = RCall.unsafe_vec(rx)
            @inbounds for cell in 1:nc, k in (Int(p[cell]) + 1):Int(p[cell + 1])
                _put!(out, cell, Int(i[k]) + 1, Float64(x[k]))
            end
        end
    elseif s.p isa Ptr{RCall.RealSxp}
        GC.@preserve s begin
            v = RCall.unsafe_vec(s)                  # column-major features × cells (or a bare numeric)
            @inbounds for cell in 1:nc, f in 1:nf
                _put!(out, cell, f, v[(cell - 1) * nf + f])
            end
        end
    else                                             # integer / logical matrix: convert while copying
        M = _dense_from_r(s, nf, nc, Float64)
        @inbounds for cell in 1:nc, f in 1:nf
            _put!(out, cell, f, M[f, cell])
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# Seurat

struct RSeuratModality
    key::String                  # SCUBA key: Seurat Key() without trailing "_"
    assay::String                # assay name
    features::Vector{String}
    layers::Vector{String}
    native::Dict{String,BPCellsMatrix}   # layer => native BPCells reader
end

struct RSeuratReduction
    key::String     # Seurat Key() without trailing "_" (e.g. "UMAP")
    name::String    # reduction name (e.g. "umap")
    ndims::Int
end

"""
    RSeuratSource(obj::RObject)

A live Seurat object (v3 `Assay` or v5 `Assay5`) as a SCUBA source. Modality keys are Seurat's
assay keys without the trailing underscore (`rna`, `ab`); assay names are accepted as aliases.
Reduction keys are the reduction keys (`UMAP`, `PC`); reduction names (`umap`, `pca`) are aliases.
"""
mutable struct RSeuratSource <: AbstractSingleCellSource
    obj::RObject
    cells::Vector{String}
    modalities::OrderedDict{String,RSeuratModality}
    assay_alias::Dict{String,String}
    reductions::OrderedDict{String,RSeuratReduction}
    red_alias::Dict{String,String}
    meta::Vector{String}
    cache::SourceCache
end

_strip_key(k::AbstractString) = endswith(k, "_") ? String(k[1:end-1]) : String(k)

function RSeuratSource(obj::RObject)
    rcopy(Bool, rcall(:inherits, obj, "Seurat")) || throw(ArgumentError("object is not a Seurat object"))
    cells = _rstrings(rcall(:colnames, obj))
    default_assay = rcopy(String, R"SeuratObject::DefaultAssay($obj)")
    assays = _rstrings(R"SeuratObject::Assays($obj)")
    mods = OrderedDict{String,RSeuratModality}()
    alias = Dict{String,String}()
    for a in vcat(default_assay, filter(!=(default_assay), assays))
        assay = R"$obj[[$a]]"
        key = _strip_key(rcopy(String, R"SeuratObject::Key($assay)"))
        isempty(key) && (key = a)
        feats = _rstrings(rcall(:rownames, assay))
        layers = _rstrings(R"scuba_seurat_layers($assay)")
        native = Dict{String,BPCellsMatrix}()
        for l in layers
            m = R"scuba_seurat_layer_data($assay, $l)"
            if rcopy(Bool, rcall(:inherits, m, "IterableMatrix"))
                info = R"scuba_bp_dir($m)"
                if !rcopy(Bool, rcall(Symbol("is.null"), info))
                    d = rcopy(String, R"$info$dir"); tr = rcopy(Bool, R"$info$transpose")
                    if isdir(d)
                        # Orientation of the files relative to the layer: decided by the layer's
                        # dimensions (features × cells); objects built as t(write_matrix_dir(t(m)))
                        # store cells × features. For square layers the R transpose flag decides
                        # (files are transposed when flag ⊻ row-major storage).
                        want = (length(feats), length(cells))
                        bm = BPCellsMatrix(d)
                        swapped = (bm.nrow, bm.ncol) != want && (bm.ncol, bm.nrow) == want
                        square = bm.nrow == bm.ncol && (bm.nrow, bm.ncol) == want
                        if swapped || (square && (tr ⊻ bm.transposed))
                            bm = BPCellsMatrix(d; transpose = true)
                        end
                        if (bm.nrow, bm.ncol) == want
                            native[l] = bm
                        else
                            @warn "SCUBA: BPCells directory shape does not match the Seurat layer; reading layer through R" assay = a layer = l dir = d disk = (bm.nrow, bm.ncol) layer_shape = want
                        end
                    end
                end
            end
        end
        mods[key] = RSeuratModality(key, a, feats, layers, native)
        alias[a] = key
    end
    reds = OrderedDict{String,RSeuratReduction}()
    ralias = Dict{String,String}()
    for r in _rstrings(R"SeuratObject::Reductions($obj)")
        red = R"$obj[[$r]]"
        key = _strip_key(rcopy(String, R"SeuratObject::Key($red)"))
        nd = rcopy(Int, R"ncol(SeuratObject::Embeddings($red))")
        reds[key] = RSeuratReduction(key, r, nd)
        ralias[r] = key
    end
    meta = _rstrings(R"colnames($obj[[]])")
    return RSeuratSource(obj, cells, mods, alias, reds, ralias, meta, SourceCache())
end

source_name(::RSeuratSource) = "RSeuratSource"
cell_ids(s::RSeuratSource) = s.cells
modality_keys(s::RSeuratSource) = collect(keys(s.modalities))
default_layer(::RSeuratSource) = "data"
canonical_modality(s::RSeuratSource, name::AbstractString) = get(s.assay_alias, String(name), String(name))
canonical_reduction(s::RSeuratSource, name::AbstractString) = get(s.red_alias, String(name), String(name))
function _mod(s::RSeuratSource, m)
    k = canonical_modality(s, m)
    haskey(s.modalities, k) || throw(ModalityNotFound(String(m), modality_keys(s)))
    return s.modalities[k]
end
feature_names(s::RSeuratSource, m::AbstractString) = _mod(s, m).features
layer_names(s::RSeuratSource, m::AbstractString) = _mod(s, m).layers
function _layer_obj(s::RSeuratSource, mod::RSeuratModality, layer)
    layer in mod.layers || throw(LayerNotFound(String(layer), mod.key, mod.layers))
    a = mod.assay
    return R"scuba_seurat_layer_data($(s.obj)[[$a]], $layer)"
end
function expression(s::RSeuratSource, m::AbstractString, layer::AbstractString, fidx::AbstractVector{Int}, cidx::AbstractVector{Int})
    mod = _mod(s, m)
    nat = get(mod.native, String(layer), nothing)
    nat === nothing || return pull(nat, fidx, cidx)
    mat = _layer_obj(s, mod, layer)
    f = _rint(collect(fidx)); c = _rint(collect(cidx))
    return _block_from_r(R"scuba_slice($mat, $f, $c)", length(fidx), length(cidx))
end
function expression_columns(s::RSeuratSource, m::AbstractString, layer::AbstractString, fidx::AbstractVector{Int},
                            cidx::AbstractVector{Int}, ::Type{T}) where {T<:Real}
    mod = _mod(s, m)
    nat = get(mod.native, String(layer), nothing)
    nat === nothing || return pull_columns(nat, fidx, cidx, T)
    mat = _layer_obj(s, mod, layer)
    f = _rint(collect(fidx)); c = _rint(collect(cidx))
    return _columns_from_r(R"scuba_slice($mat, $f, $c)", length(fidx), length(cidx), T)
end
function expression!(out, s::RSeuratSource, m::AbstractString, layer::AbstractString, fidx::AbstractVector{Int}, cidx::AbstractVector{Int})
    mod = _mod(s, m)
    nat = get(mod.native, String(layer), nothing)
    nat === nothing || return pull!(out, nat, fidx, cidx)
    mat = _layer_obj(s, mod, layer)
    f = _rint(collect(fidx)); c = _rint(collect(cidx))
    return _scatter_from_r!(out, R"scuba_slice($mat, $f, $c)", length(fidx), length(cidx))
end
supports_sparse(::RSeuratSource, ::AbstractString, ::AbstractString) = true
function expression_sparse(s::RSeuratSource, m::AbstractString, layer::AbstractString, fidx::AbstractVector{Int}, cidx::AbstractVector{Int})
    mod = _mod(s, m)
    nat = get(mod.native, String(layer), nothing)
    nat === nothing || return pull_sparse(nat, fidx, cidx)
    mat = _layer_obj(s, mod, layer)
    f = _rint(collect(fidx)); c = _rint(collect(cidx))
    return _sparse_block_from_r(R"scuba_slice($mat, $f, $c)", length(fidx), length(cidx))
end
metadata_names(s::RSeuratSource) = s.meta
function metadata_column(s::RSeuratSource, name::AbstractString, cidx::AbstractVector{Int})
    c = _rint(collect(cidx))
    return rcopy(R"$(s.obj)[[$name]][$c, 1]")
end
function metadata_table(s::RSeuratSource, cidx::AbstractVector{Int})
    c = _rint(collect(cidx))
    df = rcopy(DataFrame, R"$(s.obj)[[]][$c, , drop = FALSE]"; normalizenames = false)
    return df
end
reduction_keys(s::RSeuratSource) = collect(keys(s.reductions))
function _red(s::RSeuratSource, key)
    k = canonical_reduction(s, key)
    haskey(s.reductions, k) || throw(ReductionNotFound(String(key), reduction_keys(s)))
    return s.reductions[k]
end
reduction_ndims(s::RSeuratSource, key::AbstractString) = _red(s, key).ndims
reduction_search_names(s::RSeuratSource, key::AbstractString) = (r = _red(s, key); String[r.key, r.name])
function reduction(s::RSeuratSource, key::AbstractString, dims::AbstractVector{Int}, cidx::AbstractVector{Int})
    r = _red(s, key).name
    c = _rint(collect(cidx)); d = _rint(collect(dims))
    return Matrix{Float64}(rcopy(R"SeuratObject::Embeddings($(s.obj)[[$r]])[$c, $d, drop = FALSE]"))
end
set_metadata!(::RSeuratSource, table) =
    throw(ReadOnlySource("RSeuratSource (assign metadata in R: object@meta.data <- table)"))

# ---------------------------------------------------------------------------
# SingleCellExperiment

struct RSCEModality
    key::String
    alt::Union{Nothing,String}      # nothing = main experiment, else altExp name
    features::Vector{String}
    layers::Vector{String}
end

"""
    RSCESource(sce::RObject)

A live SingleCellExperiment (in memory or HDF5Array-backed) as a SCUBA source. Modality keys are
`mainExpName` and `altExpNames`; reduction keys are `reducedDimNames`. The default layer is
`logcounts` when present.
"""
mutable struct RSCESource <: AbstractSingleCellSource
    obj::RObject
    cells::Vector{String}
    modalities::OrderedDict{String,RSCEModality}
    reductions::OrderedDict{String,Int}
    meta::Vector{String}
    default_layer::String
    cache::SourceCache
end

function RSCESource(obj::RObject)
    rcopy(Bool, rcall(:inherits, obj, "SingleCellExperiment")) || throw(ArgumentError("object is not a SingleCellExperiment"))
    cells = _rstrings(rcall(:colnames, obj))
    main = R"SingleCellExperiment::mainExpName($obj)"
    main_key = rcopy(Bool, rcall(Symbol("is.null"), main)) ? "main" : rcopy(String, main)
    mods = OrderedDict{String,RSCEModality}()
    main_layers = _rstrings(R"SummarizedExperiment::assayNames($obj)")
    mods[main_key] = RSCEModality(main_key, nothing, _rstrings(rcall(:rownames, obj)), main_layers)
    for a in _rstrings(R"SingleCellExperiment::altExpNames($obj)")
        alt = R"SingleCellExperiment::altExp($obj, $a)"
        mods[a] = RSCEModality(a, a, _rstrings(rcall(:rownames, alt)), _rstrings(R"SummarizedExperiment::assayNames($alt)"))
    end
    reds = OrderedDict{String,Int}()
    for r in _rstrings(R"SingleCellExperiment::reducedDimNames($obj)")
        reds[r] = rcopy(Int, R"ncol(SingleCellExperiment::reducedDim($obj, $r))")
    end
    meta = _rstrings(R"colnames(SummarizedExperiment::colData($obj))")
    dl = "logcounts" in main_layers ? "logcounts" : (isempty(main_layers) ? "logcounts" : main_layers[1])
    return RSCESource(obj, cells, mods, reds, meta, dl, SourceCache())
end

source_name(::RSCESource) = "RSCESource"
cell_ids(s::RSCESource) = s.cells
modality_keys(s::RSCESource) = collect(keys(s.modalities))
default_layer(s::RSCESource) = s.default_layer
function _mod(s::RSCESource, m)
    haskey(s.modalities, m) || throw(ModalityNotFound(String(m), modality_keys(s)))
    return s.modalities[m]
end
feature_names(s::RSCESource, m::AbstractString) = _mod(s, m).features
layer_names(s::RSCESource, m::AbstractString) = _mod(s, m).layers
function _assay_obj(s::RSCESource, mod::RSCEModality, layer)
    layer in mod.layers || throw(LayerNotFound(String(layer), mod.key, mod.layers))
    if mod.alt === nothing
        return R"SummarizedExperiment::assay($(s.obj), $layer)"
    else
        a = mod.alt
        return R"SummarizedExperiment::assay(SingleCellExperiment::altExp($(s.obj), $a), $layer)"
    end
end
function expression(s::RSCESource, m::AbstractString, layer::AbstractString, fidx::AbstractVector{Int}, cidx::AbstractVector{Int})
    mat = _assay_obj(s, _mod(s, m), layer)
    f = _rint(collect(fidx)); c = _rint(collect(cidx))
    return _block_from_r(R"scuba_slice($mat, $f, $c)", length(fidx), length(cidx))
end
function expression_columns(s::RSCESource, m::AbstractString, layer::AbstractString, fidx::AbstractVector{Int},
                            cidx::AbstractVector{Int}, ::Type{T}) where {T<:Real}
    mat = _assay_obj(s, _mod(s, m), layer)
    f = _rint(collect(fidx)); c = _rint(collect(cidx))
    return _columns_from_r(R"scuba_slice($mat, $f, $c)", length(fidx), length(cidx), T)
end
function expression!(out, s::RSCESource, m::AbstractString, layer::AbstractString, fidx::AbstractVector{Int}, cidx::AbstractVector{Int})
    mat = _assay_obj(s, _mod(s, m), layer)
    f = _rint(collect(fidx)); c = _rint(collect(cidx))
    return _scatter_from_r!(out, R"scuba_slice($mat, $f, $c)", length(fidx), length(cidx))
end
supports_sparse(::RSCESource, ::AbstractString, ::AbstractString) = true
function expression_sparse(s::RSCESource, m::AbstractString, layer::AbstractString, fidx::AbstractVector{Int}, cidx::AbstractVector{Int})
    mat = _assay_obj(s, _mod(s, m), layer)
    f = _rint(collect(fidx)); c = _rint(collect(cidx))
    return _sparse_block_from_r(R"scuba_slice($mat, $f, $c)", length(fidx), length(cidx))
end
metadata_names(s::RSCESource) = s.meta
function metadata_column(s::RSCESource, name::AbstractString, cidx::AbstractVector{Int})
    c = _rint(collect(cidx))
    return rcopy(R"SummarizedExperiment::colData($(s.obj))[[$name]][$c]")
end
function metadata_table(s::RSCESource, cidx::AbstractVector{Int})
    c = _rint(collect(cidx))
    return rcopy(DataFrame, R"as.data.frame(SummarizedExperiment::colData($(s.obj)))[$c, , drop = FALSE]"; normalizenames = false)
end
reduction_keys(s::RSCESource) = collect(keys(s.reductions))
function reduction_ndims(s::RSCESource, key::AbstractString)
    haskey(s.reductions, key) || throw(ReductionNotFound(String(key), reduction_keys(s)))
    return s.reductions[key]
end
function reduction(s::RSCESource, key::AbstractString, dims::AbstractVector{Int}, cidx::AbstractVector{Int})
    reduction_ndims(s, key)
    c = _rint(collect(cidx)); d = _rint(collect(dims))
    return Matrix{Float64}(rcopy(R"as.matrix(SingleCellExperiment::reducedDim($(s.obj), $key)[$c, $d, drop = FALSE])"))
end
set_metadata!(::RSCESource, table) =
    throw(ReadOnlySource("RSCESource (assign metadata in R: colData(sce) <- table)"))

# ---------------------------------------------------------------------------
# Dispatching constructor

"""
    SCUBA.rsource(obj::RObject) -> RSeuratSource | RSCESource
"""
function SCUBA.rsource(obj::RObject)
    rcopy(Bool, rcall(:inherits, obj, "Seurat")) && return RSeuratSource(obj)
    rcopy(Bool, rcall(:inherits, obj, "SingleCellExperiment")) && return RSCESource(obj)
    throw(ArgumentError("unsupported R object of class $(rcopy(String, R"class($obj)[1]"))"))
end
SCUBA.rsource(sym::Symbol) = SCUBA.rsource(reval(string(sym)))

# ---------------------------------------------------------------------------
# R-facing glue (see R/engine.R in the SCUBA R package)

module ScubaR
# JuliaCall entry points used by the R package (R/engine.R calls Main.ScubaR.<fn>). Living here
# rather than in an R string means the glue is precompiled with the extension.
using SCUBA, RCall, DataFrames, CategoricalArrays
using OrderedCollections: OrderedDict

const S = SCUBA

# Cache of wrapped sources keyed by the R object pointer. Holding the RObject keeps the SEXP
# alive (RCall preserves it), so an address cannot be reused by a different object while cached.
# R copies on modify, so a modified object arrives with a new address and gets a fresh source.
# Least-recently-used: a Shiny app serving several datasets switches between them, and
# update_object_metadata() produces new addresses, so the oldest entry is evicted rather than
# the whole cache. Wrapping costs ~0.1 s at 1M cells (cell names are read once).
const SRC_CACHE = OrderedDict{Ptr{Cvoid}, Tuple{Any, Any}}()
const SRC_CACHE_MAX = Ref(16)
set_cache_max!(n::Integer) = (SRC_CACHE_MAX[] = max(1, Int(n)); nothing)
# File-backed sources (the R package's `scuba_source` class hands over a path, not an object): a
# second LRU keyed by the normalised path. The R side keeps the metadata, so the source is wrapped
# in an empty MetadataOverlay and Julia never resolves metadata names for it.
const PATH_CACHE = OrderedDict{String, Any}()
function getsrc(path::AbstractString)
    key = String(path)
    hit = get(PATH_CACHE, key, nothing)
    if hit !== nothing
        delete!(PATH_CACHE, key); PATH_CACHE[key] = hit
        return hit
    end
    while length(PATH_CACHE) >= SRC_CACHE_MAX[]
        k = first(keys(PATH_CACHE)); old = PATH_CACHE[k]; delete!(PATH_CACHE, k)
        try close(old) catch end
    end
    src = S.MetadataOverlay(S.open_source(key); obs = DataFrame())
    PATH_CACHE[key] = src
    return src
end
path_cache_size() = length(PATH_CACHE)
function source_close(path::AbstractString)
    src = pop!(PATH_CACHE, String(path), nothing)
    src === nothing || (try close(src) catch end)
    return nothing
end
_uns_strings(src, key) = try
    inner = src isa S.MetadataOverlay ? src.inner : src
    if (inner isa S.H5ADSource || inner isa S.MuDataSource) && S.store_haskey(inner.store, "uns/" * key)
        collect(String, string.(S.read_element(inner.store, "uns/" * key)))
    else
        String[]
    end
catch
    String[]
end
# Schema and full metadata table of a file, for the R package's file-backed object.
function source_open(path::AbstractString)
    src = getsrc(path)
    inner = src.inner
    mods = collect(String, S.modality_keys(inner))
    # Dicts cross JuliaCall as named lists; nested vectors do not
    feats = Dict{String,Vector{String}}(m => collect(String, S.feature_names(inner, m)) for m in mods)
    layers = Dict{String,Vector{String}}(m => collect(String, S.layer_names(inner, m)) for m in mods)
    reds = collect(String, S.reduction_keys(inner))
    ndims = Int[S.reduction_ndims(inner, r) for r in reds]
    dred = try S.default_reduction(inner) catch; "" end
    obs = totable(S.fetch_metadata(inner))
    return (modalities = mods, features = feats, layers = layers, reductions = reds, red_ndims = ndims,
            default_reduction = dred, default_layer = String(S.default_layer(inner)),
            uns_assays = _uns_strings(inner, "cellDIVER_assays"), uns_reductions = _uns_strings(inner, "cellDIVER_reductions"),
            cells = collect(String, S.all_cells(inner)), obs = obs)
end

function getsrc(obj::RObject)
    key = Ptr{Cvoid}(obj.p)
    hit = get(SRC_CACHE, key, nothing)
    if hit !== nothing
        delete!(SRC_CACHE, key); SRC_CACHE[key] = hit     # move to most-recently-used
        return hit[2]
    end
    while length(SRC_CACHE) >= SRC_CACHE_MAX[]
        delete!(SRC_CACHE, first(keys(SRC_CACHE)))
    end
    src = S.rsource(obj)
    SRC_CACHE[key] = (obj, src)
    return src
end
function clear_cache()
    empty!(SRC_CACHE)
    for src in values(PATH_CACHE)
        try close(src) catch end
    end
    empty!(PATH_CACHE)
    return nothing
end
cache_size() = length(SRC_CACHE)

function run(f, obj, args...; kwargs...)
    src = getsrc(obj)
    result, ws = S.with_collected_warnings(() -> f(src, args...; kwargs...))
    return result, String[w.message for w in ws]
end

asvec(x::AbstractString) = String[x]
asvec(x) = x

# Hand a table to R column by column so factors, NA and types survive the JuliaCall conversion.
# Zero-copy result columns: R numeric vectors allocated here and filled in place through SCUBA's
# `expression!` sink protocol. `totable` hands the SEXPs to R (JuliaCall passes an RObject by
# pointer), so a wide result exists once, in R, instead of once in Julia and once in R.
struct RColumns
    robjs::IdDict{Vector{Float64},RObject}
end
RColumns() = RColumns(IdDict{Vector{Float64},RObject}())
function (rc::RColumns)(::Type, ncells::Int, nfeat::Int)
    cols = Vector{Vector{Float64}}(undef, nfeat)
    for j in 1:nfeat
        r = RObject(RCall.allocArray(RCall.RealSxp, ncells))
        v = RCall.unsafe_vec(r)
        fill!(v, 0.0)                    # R does not zero fresh memory; sparse readers write nonzeros only
        rc.robjs[v] = r
        cols[j] = v
    end
    return cols
end

# `cells = false` omits the cell-name vector: R already holds the names (object column names, the
# caller's `cells`, or the scuba_source metadata) and shipping 1M strings costs ~100 ms per call.
# `robjs` maps columns that already live in R memory (see `RColumns`) to their RObject.
function totable(df::DataFrame; cells::Bool = true, robjs = nothing)
    cc = S.cell_column(df)
    names_ = String[n for n in names(df) if n != cc]
    columns = Any[]
    cat_idx = Int[]; cat_codes = Any[]; cat_levels = Any[]; cat_ordered = Bool[]
    for (k, n) in enumerate(names_)
        v = df[!, n]
        if v isa CategoricalVector
            push!(cat_idx, k)
            push!(cat_codes, Int[ismissing(x) ? -1 : levelcode(x) for x in v])
            push!(cat_levels, string.(levels(v)))
            push!(cat_ordered, isordered(v))
            push!(columns, Int[0])
        elseif robjs !== nothing && v isa Vector{Float64} && haskey(robjs, v)
            push!(columns, robjs[v])      # already an R vector: passed by pointer, not copied
        else
            push!(columns, v)
        end
    end
    # the requested spelling of each column (Julia may rename, e.g. an alt-modality bare name), so R
    # can restore the caller's column order after appending its own metadata columns
    inputs = String[(cm = colmetadata(df, n, S.COLUMN_META_KEY, nothing); cm === nothing ? n : String(cm.input)) for n in names_]
    (names = names_, columns = columns, cat_idx = cat_idx, cat_codes = cat_codes, cat_levels = cat_levels,
     cat_ordered = cat_ordered, cells = cells ? String.(S.cells(df)) : nothing, inputs = inputs)
end

_eltype(s::AbstractString) = s == "Float32" ? Float32 : s == "Float64" ? Float64 :
    throw(ArgumentError("eltype must be \"Float32\" or \"Float64\", got \"$s\""))
# `allow_empty`: the R side will append metadata columns itself, so a request whose expression
# part resolves to nothing returns an empty table plus the not-found message instead of erroring.
# `eltype` is validated but no longer changes the result: expression columns are written straight
# into R double vectors (`RColumns`), so there is no Julia-side transient left to shrink.
function fetch_data(obj, vars, layer, cells, eltype = "Float64", allow_empty = false, return_cells = true)
    _eltype(eltype)
    rc = RColumns()
    try
        df, w = run(S.fetch_data, obj, asvec(vars); layer = layer, cells = cells, alloc = rc)
        return (totable(df; cells = return_cells, robjs = rc.robjs)..., warnings = w)
    catch e
        (allow_empty && e isa S.VarsNotFound) || rethrow()
        src = getsrc(obj)
        ids = cells === nothing ? collect(String, S.all_cells(src)) : collect(String, asvec(cells))
        return (totable(S.assemble_table(ids, Pair{String,Vector{Float64}}[]); cells = return_cells)..., warnings = [sprint(showerror, e)])
    end
end
function fetch_feature(obj, features, layer, cells, modality = nothing, eltype = "Float64", return_cells = true)
    _eltype(eltype)
    rc = RColumns()
    df, w = run(S.fetch_feature, obj, asvec(features); modality = modality, layer = layer, cells = cells, alloc = rc)
    (totable(df; cells = return_cells, robjs = rc.robjs)..., warnings = w)
end
function fetch_metadata(obj, vars, cells, full_table)
    df, w = if full_table
        run(S.fetch_metadata, obj; cells = cells)
    else
        run(S.fetch_metadata, obj, asvec(vars); cells = cells)
    end
    (totable(df)..., warnings = w)
end
function fetch_reduction(obj, reduction, cells, dims, return_cells = true)
    df, w = run(S.fetch_reduction, obj, reduction; cells = cells, dims = dims)
    (totable(df; cells = return_cells)..., warnings = w)
end
# fetch_matrix for R: features × cells (the R convention), as dgCMatrix pieces (0-based i/p, x) or a
# dense matrix, plus dimnames.
function fetch_matrix(obj, features, modality, layer, cells, sparse::Bool, eltype = "Float64", return_cells = true)
    (X, cs, fs), w = run(S.fetch_matrix, obj, features === nothing ? nothing : asvec(features);
                         modality = modality, layer = layer, cells = cells, sparse = sparse, eltype = _eltype(eltype))
    if sparse
        Xt = permutedims(X)                                      # features × cells
        return (sparse = true, i = Int32.(S.SparseArrays.rowvals(Xt) .- 1), p = Int32.(Xt.colptr .- 1),
                x = S.SparseArrays.nonzeros(Xt), dim = Int[size(Xt, 1), size(Xt, 2)],
                features = collect(String, fs), cells = return_cells ? collect(String, cs) : nothing, warnings = w)
    else
        return (sparse = false, x = permutedims(X), features = collect(String, fs),
                cells = return_cells ? collect(String, cs) : nothing, warnings = w)
    end
end
fetch_cells(obj, meta_var, meta_levels) = run(S.fetch_cells, obj, meta_var, asvec(meta_levels))
get_all_cells(obj) = run(S.all_cells, obj)
meta_varnames(obj) = run(S.meta_varnames, obj)
function all_keys(obj)
    # R SCUBA returns names = assay / reduction names and values = keys; mirror that for Seurat
    src = getsrc(obj)
    ks = S.all_keys(src)
    names_ = collect(keys(ks)); vals = collect(values(ks))
    if hasproperty(src, :assay_alias)
        inv = Dict(v => k for (k, v) in src.assay_alias)
        merge!(inv, Dict(v => k for (k, v) in src.red_alias))
        names_ = [get(inv, k, k) for k in names_]
    end
    return names_, vals, String[]
end
features_in_assay(obj, assay) = run(S.features_in_modality, obj, assay)
reduction_dimnames(obj, reduction, dims) = run(S.reduction_dimnames, obj, reduction, dims)
default_layer(obj) = run(S.default_layer, obj)
default_reduction(obj) = run(S.default_reduction, obj)
function default_reduction_name(obj)
    # R SCUBA returns the reduction *name* for Seurat ("umap"), the key for others
    src = getsrc(obj)
    key = S.default_reduction(src)
    name = key
    if hasproperty(src, :reductions) && haskey(src.reductions, key) && hasproperty(src.reductions[key], :name)
        name = src.reductions[key].name
    end
    return name, String[]
end
unique_values(obj, var) = run(S.unique_values, obj, var)
n_unique(obj, var) = run(S.n_unique, obj, var)
end # module ScubaR

# ---------------------------------------------------------------------------
# Precompile workload. The extension's methods are only ever called on RObject inputs, so the
# PrecompileTools workload in SCUBA proper cannot reach them; without this the first call from
# R pays ~1.4 s of compilation. R is already running here (RCall initialised it when the
# extension loaded), so a tiny Seurat object / SingleCellExperiment is built when the packages
# are installed and every accessor is exercised once.
using PrecompileTools: @setup_workload, @compile_workload

# Module `__init__`s do not run while precompiling, so R has to be started here explicitly.
_r_running() = try
    unsafe_load(cglobal((:R_NilValue, RCall.libR), Ptr{Cvoid})) != C_NULL
catch
    false
end

function _precompile_workload()
    _r_running() || RCall.__init__()
    reval(R_HELPERS)   # normally installed by __init__, which also does not run here
    has(pkg) = rcopy(Bool, reval("requireNamespace('$pkg', quietly = TRUE)"))
    has("Matrix") || return nothing
    reval("""
    scuba_pc_counts <- Matrix::sparseMatrix(i = c(1, 2, 3, 1, 4, 2, 5), j = c(1, 1, 2, 3, 3, 4, 4), x = c(3, 1, 2, 5, 1, 4, 2),
                                            dims = c(5, 4), dimnames = list(paste0("g", 1:5), paste0("c", 1:4)))
    scuba_pc_meta <- data.frame(group = factor(c("a", "b", "a", "b")), score = c(0.5, 1.5, 2.5, 3.5), flag = c(TRUE, FALSE, TRUE, FALSE),
                                row.names = paste0("c", 1:4), stringsAsFactors = FALSE)
    scuba_pc_umap <- matrix(seq_len(8), 4, 2, dimnames = list(paste0("c", 1:4), c("UMAP_1", "UMAP_2")))
    """)
    objs = Any[]
    if has("SeuratObject")
        reval("""
        scuba_pc_seurat <- SeuratObject::CreateSeuratObject(counts = scuba_pc_counts, meta.data = scuba_pc_meta)
        scuba_pc_seurat <- SeuratObject::SetAssayData(scuba_pc_seurat, layer = "data", new.data = scuba_pc_counts)
        scuba_pc_seurat[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = scuba_pc_umap, key = "UMAP_", assay = "RNA")
        """)
        push!(objs, (reval("scuba_pc_seurat"), "rna_g1", "umap"))
        if has("BPCells")
            # cellDIVER-style objects: BPCells matrices on disk, in both storage orders
            reval("""
            scuba_pc_tmp <- tempfile("scuba_pc_bp")
            scuba_pc_bp_col <- BPCells::write_matrix_dir(methods::as(scuba_pc_counts, "IterableMatrix"), file.path(scuba_pc_tmp, "col"))
            scuba_pc_bp_row <- BPCells::transpose_storage_order(scuba_pc_bp_col, outdir = file.path(scuba_pc_tmp, "row"), tmpdir = scuba_pc_tmp)
            scuba_pc_bp_t <- BPCells::t(BPCells::write_matrix_dir(BPCells::t(methods::as(scuba_pc_counts, "IterableMatrix")), file.path(scuba_pc_tmp, "cells_x_genes")))
            scuba_pc_seurat_bp <- SeuratObject::CreateSeuratObject(counts = scuba_pc_bp_row, meta.data = scuba_pc_meta)
            scuba_pc_seurat_bp[["RNA"]]\$data <- scuba_pc_bp_col
            scuba_pc_seurat_bp[["RNA"]]\$scale.data <- scuba_pc_bp_t
            scuba_pc_seurat_bp[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = scuba_pc_umap, key = "UMAP_", assay = "RNA")
            """)
            push!(objs, (reval("scuba_pc_seurat_bp"), "rna_g1", "umap"))
        end
    end
    if has("SingleCellExperiment")
        reval("""
        scuba_pc_sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = scuba_pc_counts, logcounts = scuba_pc_counts),
                                                                   colData = scuba_pc_meta)
        SingleCellExperiment::reducedDim(scuba_pc_sce, "UMAP") <- scuba_pc_umap
        SingleCellExperiment::mainExpName(scuba_pc_sce) <- "RNA"
        """)
        push!(objs, (reval("scuba_pc_sce"), "RNA_g1", "UMAP"))
    end
    for (obj, feat, red) in objs
        src = SCUBA.rsource(obj)
        SCUBA.get_all_cells(src)
        SCUBA.meta_varnames(src)
        SCUBA.all_keys(src)
        SCUBA.features_in_assay(src)
        SCUBA.default_layer(src); SCUBA.default_reduction(src)
        SCUBA.fetch_metadata(src)
        SCUBA.fetch_metadata(src, ["group", "score"]; cells = ["c2", "c1"])
        SCUBA.fetch_metadata(src, "flag")
        SCUBA.unique_values(src, "group"); SCUBA.n_unique(src, "group")
        SCUBA.fetch_cells(src, "group", ["a"])
        SCUBA.fetch_reduction(src, red; dims = [1, 2])
        SCUBA.fetch_reduction(src, red; dims = [1, 2], cells = ["c3", "c1"])
        SCUBA.reduction_dimnames(src, red, [1, 2])
        SCUBA.fetch_data(src, [feat])
        SCUBA.fetch_data(src, [feat, "group", "UMAP_1"]; cells = ["c4", "c2"])
        SCUBA.fetch_data(src, [feat]; layer = "counts")
        SCUBA.fetch_data(src, [feat]; eltype = Float32)
        SCUBA.fetch_feature(src, ["g1", "g2"])
        SCUBA.fetch_matrix(src, ["g1", "g2"]; sparse = true)
        SCUBA.with_collected_warnings(() -> SCUBA.fetch_data(src, [feat, "nope"]))
        # the R package's entry points (objects arrive from JuliaCall exactly as from reval)
        ScubaR.fetch_data(obj, [feat, "group"], nothing, nothing)
        ScubaR.fetch_data(obj, feat, "counts", ["c2", "c1"])
        ScubaR.fetch_feature(obj, ["g1", "g2"], nothing, nothing)
        ScubaR.fetch_metadata(obj, nothing, nothing, true)
        ScubaR.fetch_metadata(obj, ["group", "score", "flag"], ["c1", "c3"], false)
        ScubaR.fetch_reduction(obj, red, nothing, [1, 2])
        ScubaR.fetch_cells(obj, "group", "a")
        ScubaR.get_all_cells(obj); ScubaR.meta_varnames(obj); ScubaR.all_keys(obj)
        ScubaR.features_in_assay(obj, nothing); ScubaR.reduction_dimnames(obj, red, [1, 2])
        ScubaR.default_layer(obj); ScubaR.default_reduction(obj); ScubaR.default_reduction_name(obj)
        ScubaR.unique_values(obj, "group"); ScubaR.n_unique(obj, "group")
        ScubaR.fetch_data(obj, feat, nothing, nothing, "Float32", false)
        ScubaR.fetch_data(obj, feat, nothing, nothing, "Float64", false, false)
        ScubaR.fetch_feature(obj, ["g1"], nothing, nothing, nothing, "Float64", false)
        ScubaR.fetch_reduction(obj, red, nothing, [1, 2], false)
        ScubaR.clear_cache()
    end
    # file-backed path (the R package's scuba_source class)
    tmp = joinpath(mktempdir(), "tiny.h5ad")
    SCUBA._write_tiny_h5ad(tmp, ["c1", "c2", "c3", "c4"], ["g1", "g2"], sparse(Float64[1 0; 0 2; 3 0; 0 4]),
                           DataFrame(group = categorical(["a", "b", "a", "b"]), score = [1.0, 2.0, 3.0, 4.0]),
                           Float64[1 2; 3 4; 5 6; 7 8])
    ScubaR.source_open(tmp)
    ScubaR.fetch_data(tmp, ["X_g1", "g2"], nothing, nothing, "Float32", false)
    ScubaR.fetch_data(tmp, ["nope"], nothing, ["c2", "c1"], "Float64", true)
    ScubaR.fetch_feature(tmp, ["g1"], nothing, nothing, nothing, "Float64")
    ScubaR.fetch_reduction(tmp, "X_umap", nothing, [1, 2])
    ScubaR.fetch_data(tmp, ["X_g1"], nothing, nothing, "Float64", false, false)
    ScubaR.fetch_reduction(tmp, "X_umap", nothing, [1, 2], false)
    ScubaR.source_close(tmp)
    reval("if (exists('scuba_pc_tmp')) unlink(scuba_pc_tmp, recursive = TRUE); rm(list = ls(pattern = '^scuba_pc_'))")
    return nothing
end

@setup_workload begin
    @compile_workload begin
        try
            _precompile_workload()
        catch err
            @warn "SCUBARCallExt precompile workload skipped" exception = (err, catch_backtrace())
        end
    end
end

end # module
