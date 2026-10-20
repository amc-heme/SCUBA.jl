# In-memory container and `snapshot`.

"""
    LazyMatrix

A matrix-like layer whose values live elsewhere (disk, another process). Implement
`Base.size(m) -> (ncells, nfeatures)` and `pull(m, feature_idx, cell_idx) -> Matrix{Float64}`
(cells × features); optionally `pull_sparse`.
"""
abstract type LazyMatrix end
pull(m::LazyMatrix, feature_idx, cell_idx) = error("pull not implemented for $(typeof(m))")
pull_sparse(m::LazyMatrix, feature_idx, cell_idx) = sparse(pull(m, feature_idx, cell_idx))

const LayerData = Union{AbstractMatrix,LazyMatrix}

"""
    Modality(features, layers)

One measured modality: feature names and a dictionary of layers, each stored
cells × features (an `AbstractMatrix`, preferably `SparseMatrixCSC`, or a `LazyMatrix`).
"""
struct Modality
    features::Vector{String}
    layers::OrderedDict{String,LayerData}
    function Modality(features::AbstractVector{<:AbstractString}, layers)
        feats = collect(String, features)
        od = OrderedDict{String,LayerData}()
        for (k, v) in pairs(layers)
            od[String(k)] = v
        end
        new(feats, od)
    end
end

"""
    MemoryAnnData(; cells, obs=DataFrame(), modalities, reductions=OrderedDict(), default_layer="X")

In-memory single-cell container. `modalities` maps modality key => `Modality` (first = main).
`reductions` maps key => cells × dims matrix. `obs` is aligned with `cells` and has no cell column.
"""
mutable struct MemoryAnnData <: AbstractSingleCellSource
    cells::Vector{String}
    obs::DataFrame
    modalities::OrderedDict{String,Modality}
    reductions::OrderedDict{String,Matrix{Float64}}
    default_layer::String
    cache::SourceCache
end

function MemoryAnnData(; cells, modalities, obs = DataFrame(), reductions = OrderedDict{String,Matrix{Float64}}(),
                       default_layer = "X")
    cells_ = collect(String, cells)
    n = length(cells_)
    mods = OrderedDict{String,Modality}()
    for (k, m) in pairs(modalities)
        mods[String(k)] = m
        for (ln, L) in m.layers
            size(L) == (n, length(m.features)) ||
                throw(DimensionMismatch("layer '$ln' of modality '$k' is $(size(L)); expected ($n, $(length(m.features)))"))
        end
    end
    isempty(mods) && throw(ArgumentError("at least one modality is required"))
    reds = OrderedDict{String,Matrix{Float64}}()
    for (k, R) in pairs(reductions)
        size(R, 1) == n || throw(DimensionMismatch("reduction '$k' has $(size(R, 1)) rows; expected $n"))
        reds[String(k)] = Matrix{Float64}(R)
    end
    obs_ = DataFrame(obs; copycols = false)
    (ncol(obs_) == 0 || nrow(obs_) == n) || throw(DimensionMismatch("obs has $(nrow(obs_)) rows; expected $n"))
    return MemoryAnnData(cells_, obs_, mods, reds, String(default_layer), SourceCache())
end

"""
    MemoryAnnData(X; cells, features, obs, layers, key="X", reductions, default_layer="X")

Single-modality convenience constructor. `X` is cells × features and becomes layer `"X"`;
`layers` adds more (name => matrix).
"""
function MemoryAnnData(X::AbstractMatrix; cells, features, obs = DataFrame(), layers = (;), key = "X",
                       reductions = OrderedDict{String,Matrix{Float64}}(), default_layer = "X")
    ls = OrderedDict{String,LayerData}("X" => X)
    for (k, v) in pairs(layers)
        ls[String(k)] = v
    end
    MemoryAnnData(; cells, obs, modalities = OrderedDict(String(key) => Modality(features, ls)), reductions, default_layer)
end

source_name(::MemoryAnnData) = "MemoryAnnData"
cell_ids(src::MemoryAnnData) = src.cells
modality_keys(src::MemoryAnnData) = collect(keys(src.modalities))
default_layer(src::MemoryAnnData) = src.default_layer
function _modality(src::MemoryAnnData, m)
    haskey(src.modalities, m) || throw(ModalityNotFound(String(m), modality_keys(src)))
    return src.modalities[m]
end
feature_names(src::MemoryAnnData, m::AbstractString) = _modality(src, m).features
layer_names(src::MemoryAnnData, m::AbstractString) = collect(keys(_modality(src, m).layers))
function _layer(src::MemoryAnnData, m, layer)
    mod = _modality(src, m)
    haskey(mod.layers, layer) || throw(LayerNotFound(String(layer), String(m), collect(keys(mod.layers))))
    return mod.layers[layer]
end
function expression(src::MemoryAnnData, m::AbstractString, layer::AbstractString,
                    fidx::AbstractVector{Int}, cidx::AbstractVector{Int})
    L = _layer(src, m, layer)
    L isa LazyMatrix ? pull(L, fidx, cidx) : gather_cells_by_features(L, cidx, fidx)
end
supports_sparse(src::MemoryAnnData, m::AbstractString, layer::AbstractString) = true
function expression_sparse(src::MemoryAnnData, m::AbstractString, layer::AbstractString,
                           fidx::AbstractVector{Int}, cidx::AbstractVector{Int})
    L = _layer(src, m, layer)
    L isa LazyMatrix && return pull_sparse(L, fidx, cidx)
    L isa SparseMatrixCSC && return gather_cf_sparse(L, cidx, fidx)
    return sparse(gather_dense(L, cidx, fidx))
end
metadata_names(src::MemoryAnnData) = names(src.obs)
metadata_column(src::MemoryAnnData, name::AbstractString, cidx::AbstractVector{Int}) = src.obs[cidx, name]
metadata_table(src::MemoryAnnData, cidx::AbstractVector{Int}) = src.obs[cidx, :]
reduction_keys(src::MemoryAnnData) = collect(keys(src.reductions))
function reduction_ndims(src::MemoryAnnData, key::AbstractString)
    haskey(src.reductions, key) || throw(ReductionNotFound(String(key), reduction_keys(src)))
    return size(src.reductions[key], 2)
end
function reduction(src::MemoryAnnData, key::AbstractString, dims::AbstractVector{Int}, cidx::AbstractVector{Int})
    haskey(src.reductions, key) || throw(ReductionNotFound(String(key), reduction_keys(src)))
    return src.reductions[key][cidx, dims]
end
function set_metadata!(src::MemoryAnnData, table)
    df = DataFrame(table; copycols = true)
    cc = metadata(df, CELL_COLUMN_KEY, nothing)
    cc === nothing || select!(df, Not(cc))
    nrow(df) == length(src.cells) || throw(DimensionMismatch("table has $(nrow(df)) rows; object has $(length(src.cells)) cells"))
    src.obs = df
    invalidate!(src)
    return src
end

# ---------------------------------------------------------------------------
# snapshot

"""
    snapshot(src; layers=:default, modalities=:all, reductions=:all, chunk=2048) -> MemoryAnnData

Materialise a source into memory for repeated queries. `layers` is `:default` (the default
layer of each modality, or its only layer), `:all`, or a vector of layer names (those
present in each modality are copied). Matrices are stored sparse when the source supports
sparse extraction, dense otherwise.
"""
function snapshot(src::AbstractSingleCellSource; layers = :default, modalities = :all, reductions = :all,
                  chunk::Int = 2048)
    n = ncells(src)
    cidx = collect(1:n)
    mods = OrderedDict{String,Modality}()
    wanted_mods = modalities === :all ? modality_keys(src) : collect(String, modalities)
    for m in wanted_mods
        check_modality(src, m)
        avail = layer_names(src, m)
        sel = if layers === :all
            avail
        elseif layers === :default
            dl = default_layer(src)
            dl in avail ? [dl] : (length(avail) == 1 ? avail : String[])
        else
            [l for l in layers if l in avail]
        end
        feats = collect(String, feature_names(src, m))
        ls = OrderedDict{String,LayerData}()
        for l in sel
            ls[String(l)] = _materialize_layer(src, m, l, length(feats), cidx, chunk)
        end
        mods[String(m)] = Modality(feats, ls)
    end
    reds = OrderedDict{String,Matrix{Float64}}()
    wanted_reds = reductions === :all ? reduction_keys(src) : collect(String, reductions)
    for r in wanted_reds
        reds[String(r)] = reduction(src, r, collect(1:reduction_ndims(src, r)), cidx)
    end
    obs = metadata_table(src, cidx)
    return MemoryAnnData(; cells = all_cells(src), obs, modalities = mods, reductions = reds,
                         default_layer = default_layer(src))
end

function _materialize_layer(src, m, l, nfeat, cidx, chunk)
    if supports_sparse(src, m, l)
        blocks = SparseMatrixCSC{Float64,Int}[]
        for start in 1:chunk:nfeat
            fidx = collect(start:min(start + chunk - 1, nfeat))
            push!(blocks, SparseMatrixCSC{Float64,Int}(expression_sparse(src, m, l, fidx, cidx)))
        end
        return length(blocks) == 1 ? blocks[1] : reduce(hcat, blocks)
    else
        out = Matrix{Float64}(undef, length(cidx), nfeat)
        for start in 1:chunk:nfeat
            fidx = collect(start:min(start + chunk - 1, nfeat))
            out[:, fidx] = expression(src, m, l, fidx, cidx)
        end
        return out
    end
end
