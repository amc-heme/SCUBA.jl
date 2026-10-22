# Public accessor API, implemented once over the AbstractSingleCellSource interface.

"""Dense-result size (bytes) above which `fetch_data` emits a `:large_result` warning (or errors, see [`DENSE_LIMIT_MODE`](@ref))."""
const MAX_DENSE_BYTES = Ref{Int}(2^30)
"""`:warn` (default) or `:error`: what `fetch_data` does when the dense result exceeds `MAX_DENSE_BYTES[]`."""
const DENSE_LIMIT_MODE = Ref{Symbol}(:warn)
"""Results below this many bytes are gathered into one dense block and split into columns (fastest);
above it, columns are filled directly so no second copy of a large result is ever held."""
const COLUMN_PATH_BYTES = Ref{Int}(256 * 2^20)

"""Thrown when `DENSE_LIMIT_MODE[] == :error` and a request exceeds `MAX_DENSE_BYTES[]`."""
struct ResultTooLarge <: SCUBAError
    bytes::Int
    limit::Int
    nfeatures::Int
    ncells::Int
end
Base.showerror(io::IO, e::ResultTooLarge) =
    print(io, "Dense result of ", e.nfeatures, " features × ", e.ncells, " cells (~", round(e.bytes / 2^30; digits = 1),
          " GiB) exceeds SCUBA.MAX_DENSE_BYTES (", round(e.limit / 2^30; digits = 1), " GiB). Request fewer features, a cell subset, ",
          "eltype = Float32, or use fetch_matrix(...; sparse = true).")

# ---------------------------------------------------------------------------
# Cells

"""
    resolve_cells(src, cells) -> (cell_idx::Vector{Int}, cell_names::Vector{String})

`cells === nothing` selects every cell in object order. Unknown IDs raise `CellsNotFound`.
"""
function resolve_cells(src::AbstractSingleCellSource, cells)
    if cells === nothing
        ids = cell_ids(src)
        # the source's own id vector is reused (no copy); result tables are built with copycols = false
        return collect(1:length(ids)), ids isa Vector{String} ? ids : collect(String, ids)
    end
    cells isa AbstractString && (cells = [cells])
    idx = cell_index(src)
    out = Vector{Int}(undef, length(cells))
    missing_cells = String[]
    for (k, c) in enumerate(cells)
        i = get(idx, String(c), 0)
        i == 0 ? push!(missing_cells, String(c)) : (out[k] = i)
    end
    isempty(missing_cells) || throw(CellsNotFound(missing_cells, length(cells)))
    return out, collect(String, cells)
end

"""`all_cells(src)` — every cell identifier, in object order (R: `get_all_cells`)."""
all_cells(src::AbstractSingleCellSource) = collect(String, cell_ids(src))
const get_all_cells = all_cells

# ---------------------------------------------------------------------------
# fetch_data

"""
    fetch_data(src, vars; layer=nothing, cells=nothing, strict=false) -> DataFrame

Fetch any mix of features, metadata columns and reduction coordinates into one table with
one row per cell (first column holds the cell IDs; see [`cell_column`](@ref)) and one column
per found variable, in input order.

- Features are addressed as `key_feature` (`RNA_GAPDH`, `X_GAPDH`, `protein_CD9-AB`) or, for
  the main modality, by bare name.
- Reduction coordinates are addressed as `reduction_dim` (`UMAP_1`, `X_umap_2`).
- Metadata columns are addressed by name, or as `obs_name`.

`layer` selects the matrix within each modality (default: [`default_layer`](@ref)). `cells`
subsets and orders the rows. With `strict=false` (default) missing variables produce a
warning; with `strict=true` they raise `VarsNotFound`. `eltype` is the element type of the
expression columns (`Float64` default; `Float32` halves the memory of large results). With
`sparse=true` expression columns are `SparseVector`s (metadata and reductions stay dense) and
the dense-size guard does not apply; see also [`eachchunk`](@ref) for very wide results.
"""
function fetch_data(src::AbstractSingleCellSource, vars::AbstractVector{<:AbstractString};
                    layer = nothing, cells = nothing, strict::Bool = false, eltype::Type{<:Real} = Float64,
                    sparse::Bool = false)
    _check_eltype(src, eltype, sparse)
    res = resolve_vars(src, vars; layer)
    if strict && !isempty(res.missing)
        throw(VarsNotFound(res.missing, length(vars), case_suggestions(src, res.missing)))
    end
    cell_idx, cell_names = resolve_cells(src, cells)
    warnings = copy(res.warnings)
    sparse || _dense_guard!(warnings, eltype, count(rv -> rv.kind == :expression, res.found), length(cell_idx))
    scuba_warn(warnings)
    return _fetch_resolved(src, res.found, cell_idx, cell_names, eltype, sparse)
end

function _check_eltype(src, ::Type{T}, sparse::Bool) where {T}
    sparse && !(T <: AbstractFloat) && throw(ArgumentError("sparse results need a floating-point eltype, got $T"))
    check_eltype(src, T)
end

# Size guard for dense results: warn (default) or throw above MAX_DENSE_BYTES.
function _dense_guard!(warnings, ::Type{T}, nfeat::Int, ncells::Int) where {T}
    est = sizeof(T) * nfeat * ncells
    if est > MAX_DENSE_BYTES[]
        DENSE_LIMIT_MODE[] == :error && throw(ResultTooLarge(est, MAX_DENSE_BYTES[], nfeat, ncells))
        push!(warnings, SCUBAWarning(:large_result,
            "The requested result is large (~$(round(est / 2^30; digits = 1)) GiB dense for $nfeat features × " *
            "$ncells cells). Consider fetch_matrix(...; sparse=true), fetch_data(...; sparse=true), eachchunk, " *
            "eltype = Float32, or fewer features."))
    end
    return warnings
end

# Build the result table for already-resolved variables over `cell_idx`. `prefetched` optionally
# supplies metadata / reduction columns (by position in `found`) already restricted to `cell_idx`;
# `eachchunk` uses it so those columns are read once rather than once per chunk.
function _fetch_resolved(src::AbstractSingleCellSource, found::AbstractVector{ResolvedVar}, cell_idx::AbstractVector{Int},
                         cell_names::AbstractVector{<:AbstractString}, eltype::Type, sparse::Bool;
                         prefetched::Union{Nothing,Dict{Int,AbstractVector}} = nothing)
    columns = Vector{AbstractVector}(undef, length(found))
    expr_groups = OrderedDict{Tuple{String,String},Vector{Int}}()
    red_groups = OrderedDict{String,Vector{Int}}()
    for (p, rv) in enumerate(found)
        if prefetched !== nothing && haskey(prefetched, p)
            columns[p] = prefetched[p]
        elseif rv.kind == :expression
            push!(get!(expr_groups, (rv.target, rv.layer::String), Int[]), p)
        elseif rv.kind == :reduction
            push!(get!(red_groups, rv.target, Int[]), p)
        else
            columns[p] = metadata_column(src, rv.target, cell_idx)
        end
    end
    for ((mod, lay), ps) in expr_groups
        fidx = [found[p].index for p in ps]
        if sparse
            X = expression_sparse(src, mod, lay, fidx, cell_idx, eltype)
            size(X) == (length(cell_idx), length(fidx)) ||
                error("backend $(source_name(src)) returned a $(size(X)) sparse block for $(length(cell_idx)) cells × $(length(fidx)) features")
            cols = _sparse_columns(X)
        elseif sizeof(eltype) * length(fidx) * length(cell_idx) <= COLUMN_PATH_BYTES[]
            M = expression(src, mod, lay, fidx, cell_idx)
            size(M) == (length(cell_idx), length(fidx)) ||
                error("backend $(source_name(src)) returned a $(size(M)) block for $(length(cell_idx)) cells × $(length(fidx)) features")
            cols = Vector{eltype}[Vector{eltype}(view(M, :, j)) for j in 1:size(M, 2)]
        else
            cols = expression_columns(src, mod, lay, fidx, cell_idx, eltype)
            length(cols) == length(fidx) && all(c -> length(c) == length(cell_idx), cols) ||
                error("backend $(source_name(src)) returned $(length(cols)) columns for $(length(fidx)) features × $(length(cell_idx)) cells")
        end
        for (j, p) in enumerate(ps)
            columns[p] = cols[j]
        end
    end
    for (key, ps) in red_groups
        dims = [found[p].index for p in ps]
        M = reduction(src, key, dims, cell_idx)
        for (j, p) in enumerate(ps)
            columns[p] = M[:, j]
        end
    end
    cols = [found[p].output => columns[p] for p in eachindex(found)]
    return assemble_table(cell_names, cols, found)
end
fetch_data(src::AbstractSingleCellSource, var::AbstractString; kw...) = fetch_data(src, [var]; kw...)

# ---------------------------------------------------------------------------
# fetch_feature / fetch_matrix: single-modality fast path

# One (modality, layer, feature indices, requested names) group per modality touched by `features`.
# With `modality` given only that modality is searched (bare or `key_`-prefixed names). Without it,
# a name prefixed with any modality's key (longest key first) goes to that modality and a bare name
# to the main modality, matching the R package's `fetch_feature`.
function _resolve_feature_groups(src, features, modality, layer)
    mods = modality === nothing ? collect(String, modality_keys(src)) : [check_modality(src, modality)]
    default = modality === nothing ? default_modality(src) : mods[1]
    sort!(mods; by = ncodeunits, rev = true)
    groups = Dict{String,Tuple{Vector{Int},Vector{String}}}()
    order = String[]
    seen = Set{Tuple{String,Int}}()
    missing_f = String[]
    for f in features
        s = String(f)
        mod = ""; i = 0
        for m in mods
            if ncodeunits(s) > ncodeunits(m) + 1 && startswith(s, m) && codeunit(s, ncodeunits(m) + 1) == UInt8('_')
                i = get(feature_index(src, m), strip_key(s, m), 0)
                i == 0 || (mod = m; break)
            end
        end
        if i == 0
            i = get(feature_index(src, default), s, 0)
            mod = default
        end
        if i == 0
            push!(missing_f, s)
        elseif !((mod, i) in seen)
            push!(seen, (mod, i))
            g = get!(groups, mod) do
                push!(order, mod); (Int[], String[])
            end
            push!(g[1], i); push!(g[2], s)
        end
    end
    isempty(order) && throw(VarsNotFound(missing_f, length(features), case_suggestions(src, missing_f)))
    isempty(missing_f) || scuba_warn(:some_vars_not_found,
        "The following requested features were not found" *
        (modality === nothing ? "" : " in modality '$default'") * _shown(missing_f) * ": " *
        _join(missing_f), missing_f...)
    return [(m, resolve_layer(src, m, layer), groups[m][1], groups[m][2]) for m in order]
end

# Single-modality resolution for `fetch_matrix`.
function _resolve_features(src, features, modality, layer)
    groups = _resolve_feature_groups(src, features, modality, layer)
    length(groups) == 1 || throw(ArgumentError("features span several modalities (" *
        join((g[1] for g in groups), ", ") * "); fetch_matrix takes features from one modality"))
    return groups[1]
end

"""
    fetch_feature(src, features; modality=nothing, layer=nothing, cells=nothing) -> DataFrame

Fetch expression features. Bare names are looked up in `modality` (default: the main modality);
names prefixed with a modality's key (`ab_CD34`) come from that modality, so one call may span
several modalities. Columns keep the requested names. Unlike `fetch_data`, metadata and reductions
are not searched.
"""
function fetch_feature(src::AbstractSingleCellSource, features::AbstractVector{<:AbstractString};
                       modality = nothing, layer = nothing, cells = nothing, eltype::Type{<:Real} = Float64,
                       sparse::Bool = false)
    _check_eltype(src, eltype, sparse)
    groups = _resolve_feature_groups(src, features, modality, layer)
    cell_idx, cell_names = resolve_cells(src, cells)
    pairs = Pair{String,AbstractVector}[]
    resolved = ResolvedVar[]
    for (mod, lay, idx, names) in groups
        cols = sparse ? _sparse_columns(expression_sparse(src, mod, lay, idx, cell_idx, eltype)) :
                        expression_columns(src, mod, lay, idx, cell_idx, eltype)
        for (j, n) in enumerate(names)
            push!(pairs, n => cols[j]); push!(resolved, ResolvedVar(n, n, :expression, mod, lay, idx[j]))
        end
    end
    if length(groups) > 1   # restore the requested order across modalities
        pos = Dict(String(f) => k for (k, f) in enumerate(features))
        perm = sortperm([pos[p.first] for p in pairs])
        pairs = pairs[perm]; resolved = resolved[perm]
    end
    return assemble_table(cell_names, pairs, resolved)
end
fetch_feature(src::AbstractSingleCellSource, feature::AbstractString; kw...) = fetch_feature(src, [feature]; kw...)

"""
    fetch_matrix(src, features; modality=nothing, layer=nothing, cells=nothing, sparse=false)
        -> (X, cells, features)

Expression block (cells × features) as a plain `Matrix{eltype}` or, with `sparse=true`, a
`SparseMatrixCSC{eltype}`. Intended for statistics-heavy consumers (dot plots, heatmaps).
`features === nothing` selects every feature of the modality. Dense results are subject to the
`MAX_DENSE_BYTES` guard; sparse results are not.
"""
function fetch_matrix(src::AbstractSingleCellSource, features = nothing;
                      modality = nothing, layer = nothing, cells = nothing, sparse::Bool = false,
                      eltype::Type{<:Real} = Float64)
    _check_eltype(src, eltype, sparse)
    if features === nothing
        mod = modality === nothing ? default_modality(src) : check_modality(src, modality)
        lay = resolve_layer(src, mod, layer)
        names = collect(String, feature_names(src, mod))
        idx = collect(1:length(names))
    else
        mod, lay, idx, names = _resolve_features(src, features, modality, layer)
    end
    cell_idx, cell_names = resolve_cells(src, cells)
    if sparse
        X = expression_sparse(src, mod, lay, idx, cell_idx, eltype)
    else
        scuba_warn(_dense_guard!(SCUBAWarning[], eltype, length(idx), length(cell_idx)))
        X = _dense_block(src, mod, lay, idx, cell_idx, eltype)
    end
    return X, cell_names, names
end

# ---------------------------------------------------------------------------
# Metadata

"""`meta_varnames(src)` — names of the metadata columns."""
meta_varnames(src::AbstractSingleCellSource) = collect(String, metadata_names(src))

function _check_meta_vars(src, vars)
    avail = metadata_names(src)
    aset = Set{String}(String.(avail))
    missing_v = [String(v) for v in vars if !(String(v) in aset)]
    isempty(missing_v) || throw(MetadataNotFound(missing_v, collect(String, avail)))
    return nothing
end

"""
    fetch_metadata(src, vars; cells=nothing) -> DataFrame
    fetch_metadata(src, var; cells=nothing) -> AbstractVector
    fetch_metadata(src; cells=nothing) -> DataFrame   (full table)

Metadata columns for the requested cells. Any unknown variable raises `MetadataNotFound`.
"""
function fetch_metadata(src::AbstractSingleCellSource, vars::AbstractVector{<:AbstractString}; cells = nothing)
    _check_meta_vars(src, vars)
    cell_idx, cell_names = resolve_cells(src, cells)
    uvars = unique(String.(vars))
    cols = [v => metadata_column(src, v, cell_idx) for v in uvars]
    resolved = [ResolvedVar(v, v, :metadata, v, nothing, 0) for v in uvars]
    return assemble_table(cell_names, cols, resolved)
end
function fetch_metadata(src::AbstractSingleCellSource, var::AbstractString; cells = nothing)
    _check_meta_vars(src, [var])
    cell_idx, _ = resolve_cells(src, cells)
    return metadata_column(src, var, cell_idx)
end
function fetch_metadata(src::AbstractSingleCellSource; cells = nothing)
    cell_idx, cell_names = resolve_cells(src, cells)
    tbl = metadata_table(src, cell_idx)
    cols = [String(n) => tbl[!, n] for n in names(tbl)]
    resolved = [ResolvedVar(String(n), String(n), :metadata, String(n), nothing, 0) for n in names(tbl)]
    return assemble_table(cell_names, cols, resolved)
end

_is_categorical(col::CategoricalArray) = true
_is_categorical(col::AbstractVector) = nonmissingtype(eltype(col)) <: Union{AbstractString,Bool,Symbol,Char}

_level_string(x) = ismissing(x) ? missing : string(x)
_level_string(x::CategoricalValue) = string(unwrap(x))

"""
    fetch_cells(src, meta_var, levels) -> Vector{String}

IDs of the cells whose categorical metadata value is in `levels`, in object order.
"""
function fetch_cells(src::AbstractSingleCellSource, meta_var::AbstractString, levels)
    levels === nothing && throw(ArgumentError("`levels` must be provided"))
    _check_meta_vars(src, [meta_var])
    col = metadata_column(src, meta_var, 1:ncells(src))
    _is_categorical(col) || throw(NonCategoricalVar(String(meta_var), eltype(col)))
    wanted = Set{String}(string(l) for l in (levels isa AbstractVector ? levels : [levels]))
    ids = cell_ids(src)
    out = String[]
    for (i, x) in enumerate(col)
        s = _level_string(x)
        s === missing && continue
        s in wanted && push!(out, String(ids[i]))
    end
    return out
end

"""
    unique_values(src, var) -> AbstractVector

Unique values of a categorical metadata column, in order of first appearance. A
`CategoricalArray` column yields a `CategoricalArray` (levels preserved).
"""
function unique_values(src::AbstractSingleCellSource, var::AbstractString)
    _check_meta_vars(src, [var])
    col = metadata_column(src, var, 1:ncells(src))
    _is_categorical(col) || throw(NonCategoricalVar(String(var), eltype(col)))
    return unique(col)
end

"""`n_unique(src, var)` — number of unique values of a categorical metadata column."""
n_unique(src::AbstractSingleCellSource, var::AbstractString) = length(unique_values(src, var))

"""
    update_object_metadata!(src, table)

Replace the metadata table of `src` (alias of [`set_metadata!`](@ref)). A SCUBA cell column in
`table` is dropped; the remaining rows must be in object order.
"""
update_object_metadata!(src::AbstractSingleCellSource, table) = set_metadata!(src, table)

# ---------------------------------------------------------------------------
# Modalities and features

"""`features_in_modality(src, modality=default_modality(src))` — feature names of a modality (R: `features_in_assay`)."""
function features_in_modality(src::AbstractSingleCellSource, modality = nothing)
    mod = modality === nothing ? default_modality(src) : check_modality(src, modality)
    return collect(String, feature_names(src, mod))
end
const features_in_assay = features_in_modality

# ---------------------------------------------------------------------------
# Reductions

const DEFAULT_REDUCTION_ORDER = ("umap", "tsne", "pca")

_reduction_basename(key) = lowercase(startswith(key, "X_") ? key[3:end] : key)

"""
    default_reduction(src) -> String

The preferred reduction: the first of umap, tsne, pca present (case-insensitive, with or
without an `X_` prefix). Raises `NoDefaultReduction` if none is found.
"""
function default_reduction(src::AbstractSingleCellSource)
    reds = collect(String, reduction_keys(src))
    for want in DEFAULT_REDUCTION_ORDER
        for r in reds
            any(n -> _reduction_basename(n) == want, reduction_search_names(src, r)) && return r
        end
    end
    throw(NoDefaultReduction(reds))
end

function _check_dims(src, key, dims)
    d = collect(Int, dims)
    nd = reduction_ndims(src, key)
    all(1 .<= d .<= nd) || throw(InvalidDims(d, String(key), nd))
    return d
end

"""
    fetch_reduction(src, reduction=default_reduction(src); cells=nothing, dims=1:2) -> DataFrame

Embedding coordinates for `dims`, with the cell column first and columns named by
[`reduction_dimnames`](@ref).
"""
function fetch_reduction(src::AbstractSingleCellSource, red = nothing; cells = nothing, dims = 1:2)
    key = red === nothing ? default_reduction(src) : check_reduction(src, red)
    d = _check_dims(src, key, dims)
    cell_idx, cell_names = resolve_cells(src, cells)
    M = reduction(src, key, d, cell_idx)
    names_ = reduction_dimnames_impl(src, key, d)
    resolved = [ResolvedVar(n, n, :reduction, key, nothing, dd) for (n, dd) in zip(names_, d)]
    return assemble_table(cell_names, [n => M[:, j] for (j, n) in enumerate(names_)], resolved)
end

"""`reduction_dimnames(src, reduction, dims)` — column names used for a reduction's dimensions."""
function reduction_dimnames(src::AbstractSingleCellSource, red::AbstractString, dims)
    key = check_reduction(src, red)
    return reduction_dimnames_impl(src, key, _check_dims(src, key, dims))
end
