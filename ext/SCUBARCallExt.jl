# Live R objects (Seurat, SingleCellExperiment) as SCUBA sources, through RCall.
#
# Nothing is copied at construction except names and shapes. Expression blocks are sliced in R
# and only the requested (features × cells) rectangle crosses into Julia. BPCells layers that are
# plain on-disk matrices are read natively by `BPCellsMatrix`, bypassing R entirely.
module SCUBARCallExt

using SCUBA, RCall, DataFrames, CategoricalArrays, SparseArrays, OrderedCollections
using SCUBA: AbstractSingleCellSource, SourceCache, LazyMatrix, BPCellsMatrix, pull, pull_sparse,
             gather_features_by_cells, gather_fc_sparse, ModalityNotFound, LayerNotFound, ReductionNotFound,
             ReadOnlySource
import SCUBA: cell_ids, modality_keys, default_modality, default_layer, feature_names, layer_names,
              expression, expression_sparse, supports_sparse, metadata_names, metadata_column, metadata_table,
              reduction_keys, reduction_ndims, reduction, source_name, set_metadata!, invalidate!,
              canonical_modality, canonical_reduction, reduction_search_names

const R_HELPERS = """
scuba_bp_dir <- function(x) {
  # Return the on-disk directory of a BPCells matrix if it is a plain MatrixDir, possibly wrapped
  # only in dimension-renaming wrappers; NULL otherwise (transformed/queued matrices need R).
  while (!methods::is(x, "MatrixDir")) {
    cls <- class(x)[1]
    if (!(cls %in% c("RenameDims")) || !methods::.hasSlot(x, "matrix")) return(NULL)
    x <- x@matrix
  }
  if (!methods::.hasSlot(x, "dir")) return(NULL)
  x@dir
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
function _block_from_r(s::RObject, nf::Int, nc::Int)
    if _is_dgc(s)
        return gather_features_by_cells(_dgc_to_sparse(s), collect(1:nc), collect(1:nf))
    else
        return Matrix{Float64}(permutedims(rcopy(Matrix{Float64}, s)))
    end
end
function _sparse_block_from_r(s::RObject, nf::Int, nc::Int)
    _is_dgc(s) && return gather_fc_sparse(_dgc_to_sparse(s), collect(1:nc), collect(1:nf))
    return sparse(Matrix{Float64}(permutedims(rcopy(Matrix{Float64}, s))))
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
                dir = R"scuba_bp_dir($m)"
                if !rcopy(Bool, rcall(Symbol("is.null"), dir))
                    d = rcopy(String, dir)
                    isdir(d) && (native[l] = BPCellsMatrix(d))
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

end # module
