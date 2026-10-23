# Native reader for MuData (.h5mu) containers: one anndata group per modality under `mod/`, a
# shared `obs` table over the union of cells, and per-modality `obsmap`s (1-based row of each
# global cell in the modality, 0 when absent). Only `axis = 0` (shared observations) is supported.
#
# Naming (mirrors the SCUBA R package's MuData semantics):
#   modality key = modality name (`rna_GAPDH` reads mod/rna/X); a modality's extra anndata
#     modalities (raw, obsm data frames) are `<mod>_<key>`;
#   reductions = `<mod>_<obsm key>` (`rna_X_umap`), plus any numeric top-level `obsm` array under
#     its own name; a bare child key (`X_umap`) is an alias when exactly one modality has it and no
#     global reduction uses the name;
#   metadata = global obs columns (`mod:col` shown as `mod_col`), every modality obs column as
#     `<mod>_<col>`, and bare `<col>` for columns found in exactly one modality and not globally;
#   cells absent from a modality: NaN (dense / column results), no entry (sparse), `missing`
#     (metadata), NaN rows (reductions).

struct MetaRef
    path::String
    modality::Union{Nothing,String}
end
struct MuModality
    child::H5ADSource
    child_key::String
    obsmap::Vector{Int}          # global cell i -> child row (0 = absent)
    all_present::Bool
end
struct MuReduction
    child::Union{Nothing,H5ADSource}
    key::String
    ref::Union{Nothing,ArrayRef}             # global reduction
    obsmap::Union{Nothing,Vector{Int}}
end

"""
    MuDataSource(path) / open_h5mu(path)

Lazily read a MuData `.h5mu` file (see the module notes for how modalities, reductions and
metadata are named). Cells absent from a modality read as `NaN` / `missing`. Close with `close`.
"""
mutable struct MuDataSource <: AbstractSingleCellSource
    path::String
    store::AbstractStore
    cells::Vector{String}
    mod_order::Vector{String}
    modalities::OrderedDict{String,MuModality}
    reductions::OrderedDict{String,MuReduction}
    reduction_alias::Dict{String,String}
    meta::OrderedDict{String,MetaRef}
    cache::SourceCache
end

function MuDataSource(path::AbstractString)
    isfile(path) || throw(ArgumentError("file not found: $path"))
    file = h5open(String(path), "r")
    try
        return mudata_source_from_store(H5Store(file), String(path))
    catch
        close(file)
        rethrow()
    end
end
open_h5mu(path::AbstractString) = MuDataSource(path)
Base.close(src::MuDataSource) = (src.store isa H5Store && close(src.store.file); nothing)
Base.isopen(src::MuDataSource) = src.store isa H5Store ? isopen(src.store.file) : true

_obsmap_from_names(cells, child_cells) = begin
    pos = Dict{String,Int}(c => i for (i, c) in enumerate(child_cells))
    Int[get(pos, c, 0) for c in cells]
end

"""Build a `MuDataSource` over any `AbstractStore` holding the MuData layout."""
function mudata_source_from_store(store::AbstractStore, path::String)
    attrs = store_root_attrs(store)
    et = string(get(attrs, "encoding-type", ""))
    et == "MuData" || throw(ArgumentError("$path is not a MuData container (encoding-type '$et')"))
    axis = Int(get(attrs, "axis", 0))
    axis == 0 || throw(ArgumentError("only MuData files with axis = 0 (shared observations) are supported; this file has axis = $axis"))
    cells = dataframe_index(store, "obs")
    gcols = dataframe_columns(store, "obs")
    children = store_haskey(store, "mod") ? collect(String, store_children(store, "mod")) : String[]
    order = String.(get(store_attrs(store, "mod"), "mod-order", sort(children)))
    mod_order = vcat([m for m in order if m in children], sort(setdiff(children, order)))
    isempty(mod_order) && throw(ArgumentError("$path has no modalities under mod/"))

    modalities = OrderedDict{String,MuModality}()
    reductions = OrderedDict{String,MuReduction}()
    for m in mod_order
        child = anndata_source_from_store(store, path; prefix = "mod/" * m)
        obsmap = store_haskey(store, "obsmap/" * m) ? Int.(store_read(store, "obsmap/" * m)) :
                 _obsmap_from_names(cells, child.cells)
        length(obsmap) == length(cells) ||
            error("$path: obsmap/$m has $(length(obsmap)) entries for $(length(cells)) cells")
        allp = all(!=(0), obsmap)
        for ck in modality_keys(child)
            modalities[ck == "X" ? m : m * "_" * ck] = MuModality(child, ck, obsmap, allp)
        end
        for r in reduction_keys(child)
            reductions[m * "_" * r] = MuReduction(child, r, nothing, obsmap)
        end
    end
    if store_haskey(store, "obsm")
        for k in store_children(store, "obsm")
            k in mod_order && continue                           # boolean membership masks
            p = "obsm/" * k
            encoding_type(store, p) == "array" || continue
            sz = store_size(store, p)
            length(sz) == 2 || continue
            haskey(reductions, k) || (reductions[k] = MuReduction(nothing, k, ArrayRef(p, Int(sz[2]), Int(sz[1])), nothing))
        end
    end
    counts = Dict{String,Int}()
    for (k, r) in reductions
        r.child === nothing || (counts[r.key] = get(counts, r.key, 0) + 1)
    end
    alias = Dict{String,String}()
    for (k, r) in reductions
        r.child === nothing && continue
        counts[r.key] == 1 && !haskey(reductions, r.key) && (alias[r.key] = k)
    end

    meta = OrderedDict{String,MetaRef}()
    for c in gcols
        meta[replace(c, ":" => "_")] = MetaRef("obs/" * c, nothing)
    end
    percol = OrderedDict{String,Vector{String}}()
    for m in mod_order
        child = modalities[m].child
        for c in metadata_names(child)
            key = m * "_" * c
            haskey(meta, key) || (meta[key] = MetaRef(child.obs_path * "/" * c, m))
            push!(get!(percol, c, String[]), m)
        end
    end
    for (c, ms) in percol
        length(ms) == 1 && !haskey(meta, c) && (meta[c] = MetaRef(modalities[ms[1]].child.obs_path * "/" * c, ms[1]))
    end
    return MuDataSource(path, store, cells, mod_order, modalities, reductions, alias, meta, SourceCache())
end

source_name(::MuDataSource) = "MuDataSource"
cell_ids(src::MuDataSource) = src.cells
modality_keys(src::MuDataSource) = collect(keys(src.modalities))
default_modality(src::MuDataSource) = src.mod_order[1]
default_layer(::MuDataSource) = "X"
function _mm(src::MuDataSource, m::AbstractString)
    haskey(src.modalities, m) || throw(ModalityNotFound(String(m), modality_keys(src)))
    return src.modalities[m]
end
feature_names(src::MuDataSource, m::AbstractString) = (mm = _mm(src, m); feature_names(mm.child, mm.child_key))
layer_names(src::MuDataSource, m::AbstractString) = (mm = _mm(src, m); layer_names(mm.child, mm.child_key))
supports_sparse(src::MuDataSource, m::AbstractString, l::AbstractString) = (mm = _mm(src, m); supports_sparse(mm.child, mm.child_key, l))

check_eltype(::MuDataSource, ::Type{T}) where {T} =
    T <: AbstractFloat ? nothing : throw(ArgumentError("MuData results need a floating-point eltype (absent cells are NaN), got $T"))

# Absent cells: NaN in dense / column sinks, nothing in sparse sinks.
function _mark_absent!(out::AbstractMatrix, rows, nf::Int)
    eltype(out) <: AbstractFloat || throw(ArgumentError("MuData results need a floating-point eltype (absent cells are NaN), got $(eltype(out))"))
    @inbounds for j in 1:nf, i in rows
        out[i, j] = NaN
    end
    return out
end
function _mark_absent!(cols::Vector{<:AbstractVector}, rows, nf::Int)
    nf == 0 && return cols
    eltype(cols[1]) <: AbstractFloat || throw(ArgumentError("MuData results need a floating-point eltype (absent cells are NaN), got $(eltype(cols[1]))"))
    @inbounds for j in 1:nf, i in rows
        cols[j][i] = NaN
    end
    return cols
end
_mark_absent!(s::SparseSink, rows, nf::Int) = s
_mark_absent!(s::GroupSumSink, rows, nf::Int) = s          # absent cells add nothing to a group sum
cell_present(src::MuDataSource, m::AbstractString, cidx::AbstractVector{Int}) = _mm(src, m).obsmap[cidx] .!= 0

function expression!(out, src::MuDataSource, m::AbstractString, layer::AbstractString,
                     fidx::AbstractVector{Int}, cidx::AbstractVector{Int})
    mm = _mm(src, m)
    rows = mm.obsmap[cidx]
    if mm.all_present
        return expression!(out, mm.child, mm.child_key, layer, fidx, rows)
    end
    present = findall(!=(0), rows)
    isempty(present) || expression!(RowMapSink(out, present), mm.child, mm.child_key, layer, fidx, rows[present])
    _mark_absent!(out, findall(==(0), rows), length(fidx))
    return out
end
expression(src::MuDataSource, m::AbstractString, l::AbstractString, fidx::AbstractVector{Int}, cidx::AbstractVector{Int}) =
    expression!(zeros(Float64, length(cidx), length(fidx)), src, m, l, fidx, cidx)
expression_sparse(src::MuDataSource, m::AbstractString, l::AbstractString, fidx::AbstractVector{Int},
                  cidx::AbstractVector{Int}, ::Type{T}) where {T} =
    to_sparse_matrix(expression!(SparseSink(T, length(fidx)), src, m, l, fidx, cidx), length(cidx))
expression_sparse(src::MuDataSource, m::AbstractString, l::AbstractString, fidx::AbstractVector{Int}, cidx::AbstractVector{Int}) =
    expression_sparse(src, m, l, fidx, cidx, Float64)

metadata_names(src::MuDataSource) = collect(keys(src.meta))
function metadata_column(src::MuDataSource, name::AbstractString, cidx::AbstractVector{Int})
    haskey(src.meta, name) || throw(MetadataNotFound([String(name)], metadata_names(src)))
    ref = src.meta[name]
    col = read_element(src.store, ref.path)
    ref.modality === nothing && return col[cidx]
    rows = src.modalities[ref.modality].obsmap[cidx]
    all(!=(0), rows) && return col[rows]
    out = allowmissing(col[max.(rows, 1)])
    out[rows .== 0] .= missing
    return out
end

reduction_keys(src::MuDataSource) = collect(keys(src.reductions))
reduction_aliases(src::MuDataSource) = src.reduction_alias
canonical_reduction(src::MuDataSource, name::AbstractString) = get(src.reduction_alias, String(name), String(name))
function _mr(src::MuDataSource, key::AbstractString)
    k = canonical_reduction(src, key)
    haskey(src.reductions, k) || throw(ReductionNotFound(String(key), reduction_keys(src)))
    return k, src.reductions[k]
end
reduction_ndims(src::MuDataSource, key::AbstractString) = (r = _mr(src, key)[2]; r.child === nothing ? r.ref.k : reduction_ndims(r.child, r.key))
function reduction_search_names(src::MuDataSource, key::AbstractString)
    k, r = _mr(src, key)
    names = String[k]
    r.child === nothing || (get(src.reduction_alias, r.key, "") == k && push!(names, r.key))
    return names
end
function reduction(src::MuDataSource, key::AbstractString, dims::AbstractVector{Int}, cidx::AbstractVector{Int})
    _, r = _mr(src, key)
    r.child === nothing && return read_array_block(src.store, r.ref, dims, cidx)
    rows = r.obsmap[cidx]
    all(!=(0), rows) && return reduction(r.child, r.key, dims, rows)
    out = fill(NaN, length(cidx), length(dims))
    present = findall(!=(0), rows)
    isempty(present) || (out[present, :] = reduction(r.child, r.key, dims, rows[present]))
    return out
end
