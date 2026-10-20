# Decoders for the anndata on-disk specification (encoding-type / encoding-version attrs),
# written against the AbstractStore API.

struct MatrixRef
    kind::Symbol            # :csr | :csc | :dense | :dataframe
    path::String
    n_obs::Int
    n_var::Int
    columns::Vector{String} # for :dataframe only
end
MatrixRef(kind, path, n_obs, n_var) = MatrixRef(kind, String(path), n_obs, n_var, String[])

struct ArrayRef             # obsm numeric array: python shape (n_obs, k)
    path::String
    n_obs::Int
    k::Int
end

_attr(attrs, key, default = nothing) = get(attrs, key, default)

"""Encoding type of an element: the `encoding-type` attribute, or an inference for legacy files."""
function encoding_type(store::AbstractStore, path::AbstractString)
    attrs = store_attrs(store, path)
    et = _attr(attrs, "encoding-type")
    et === nothing || return String(et)
    if store_isgroup(store, path)
        ch = Set(store_children(store, path))
        if "codes" in ch && "categories" in ch
            return "categorical"
        elseif "data" in ch && "indices" in ch && "indptr" in ch
            fmt = _attr(attrs, "h5sparse_format")
            fmt === nothing && return "csr_matrix"
            return String(fmt) * "_matrix"
        elseif "values" in ch && "mask" in ch
            return store_eltype(store, path * "/values") <: AbstractString ? "nullable-string-array" : "nullable-integer"
        elseif "_index" in ch || haskey(attrs, "_index")
            return "dataframe"
        else
            return "dict"
        end
    else
        T = store_eltype(store, path)
        T <: AbstractString && return "string-array"
        return "array"
    end
end

# ---------------------------------------------------------------------------
# 1-D elements (data frame columns)

_to_bool(x) = x isa Bool ? x : Bool(x)

"""Read a 1-D element (a data frame column) fully."""
function read_element(store::AbstractStore, path::AbstractString)
    et = encoding_type(store, path)
    if et == "categorical"
        codes = store_read(store, path * "/codes")
        cats = store_read(store, path * "/categories")
        ordered = _to_bool(_attr(store_attrs(store, path), "ordered", false))
        return _categorical_from_codes(codes, cats, ordered)
    elseif et == "nullable-string-array"
        vals = store_read(store, path * "/values")
        mask = store_read(store, path * "/mask")
        any(_to_bool, mask) || return collect(String, string.(vals))
        out = Vector{Union{Missing,String}}(undef, length(vals))
        for i in eachindex(vals)
            out[i] = _to_bool(mask[i]) ? missing : String(string(vals[i]))
        end
        return out
    elseif et == "nullable-integer" || et == "nullable-boolean"
        vals = store_read(store, path * "/values")
        mask = store_read(store, path * "/mask")
        T = et == "nullable-boolean" ? Bool : eltype(vals)
        out = Vector{Union{Missing,T}}(undef, length(vals))
        for i in eachindex(vals)
            out[i] = _to_bool(mask[i]) ? missing : T(vals[i])
        end
        return out
    elseif et == "string-array" || et == "array" || et == "string" || et == "numeric-scalar"
        x = store_read(store, path)
        return x isa AbstractVector ? x : vec(x)
    else
        error("Unsupported anndata encoding '$et' at $path")
    end
end

function _categorical_from_codes(codes::AbstractVector{<:Integer}, cats::AbstractVector, ordered::Bool)
    T = eltype(cats) <: AbstractString ? String : eltype(cats)
    levels_ = collect(T, cats)
    anymissing = any(<(0), codes)
    V = anymissing ? Union{Missing,T} : T
    vals = Vector{V}(undef, length(codes))
    @inbounds for i in eachindex(codes)
        vals[i] = codes[i] < 0 ? missing : levels_[codes[i] + 1]
    end
    return categorical(vals; levels = levels_, ordered = ordered)
end

"""Column names of an anndata data frame element, in stored order."""
function dataframe_columns(store::AbstractStore, path::AbstractString)
    attrs = store_attrs(store, path)
    co = _attr(attrs, "column-order")
    if co !== nothing
        return collect(String, co isa AbstractString ? [co] : co)
    end
    idx = String(_attr(attrs, "_index", "_index"))
    return [c for c in store_children(store, path) if c != idx && c != "__categories"]
end

"""Index (row names) of an anndata data frame element."""
function dataframe_index(store::AbstractStore, path::AbstractString)
    attrs = store_attrs(store, path)
    idx = String(_attr(attrs, "_index", "_index"))
    col = read_element(store, path * "/" * idx)
    return String[ismissing(x) ? "" : string(x) for x in col]
end

# ---------------------------------------------------------------------------
# Matrices

"""Describe the matrix element at `path` (X, a layer, raw/X, or an obsm entry)."""
function matrix_ref(store::AbstractStore, path::AbstractString)
    et = encoding_type(store, path)
    attrs = store_attrs(store, path)
    if et == "csr_matrix" || et == "csc_matrix"
        shape = _attr(attrs, "shape", _attr(attrs, "h5sparse_shape"))
        shape === nothing && error("sparse matrix at $path has no shape attribute")
        return MatrixRef(et == "csr_matrix" ? :csr : :csc, path, Int(shape[1]), Int(shape[2]))
    elseif et == "array"
        sz = store_size(store, path)       # Julia order: (n_var, n_obs)
        length(sz) == 2 || error("array at $path is not 2-D")
        return MatrixRef(:dense, path, Int(sz[2]), Int(sz[1]))
    elseif et == "dataframe"
        cols = dataframe_columns(store, path)
        n = length(dataframe_index(store, path))
        return MatrixRef(:dataframe, path, n, length(cols), cols)
    else
        error("Element at $path (encoding '$et') is not a matrix")
    end
end

const MAX_READ_ELEMENTS = Ref{Int}(1 << 23)   # ≈ 64 MiB of Float64 per read

# Iterate over `(a, b, members)` blocks: contiguous 1-based ranges a:b of a compressed
# sparse axis covering consecutive members (sorted), split so a read stays under the limit.
function _foreach_block(f, indptr::Vector{Int}, runs::Vector{UnitRange{Int}})
    limit = MAX_READ_ELEMENTS[]
    for run in runs
        start = first(run)
        while start <= last(run)
            stop = start
            while stop < last(run) && (indptr[stop + 2] - indptr[start + 1]) <= limit
                stop += 1
            end
            a = indptr[start] + 1
            b = indptr[stop + 1]
            a <= b && f(a, b, start:stop)
            start = stop + 1
        end
    end
end

"""
    read_sparse_block(store, ref, indptr, cell_idx, feature_idx) -> Matrix{Float64}

Dense cells × features block from a CSR/CSC matrix element. `indptr` is the cached,
0-based index pointer (length n+1 along the compressed axis).
"""
function read_sparse_block(store::AbstractStore, ref::MatrixRef, indptr::Vector{Int},
                           cell_idx::AbstractVector{Int}, feature_idx::AbstractVector{Int})
    out = zeros(Float64, length(cell_idx), length(feature_idx))
    data_path = ref.path * "/data"; idx_path = ref.path * "/indices"
    if ref.kind == :csc
        # compressed axis = features (columns); indices are cell (row) ids
        pos = position_lookup(ref.n_obs, cell_idx)
        colpos = position_lookup(ref.n_var, feature_idx)
        runs = consecutive_runs(sort(feature_idx))
        _foreach_block(indptr, runs) do a, b, members
            ind = store_read(store, idx_path, a:b)
            val = store_read(store, data_path, a:b)
            for f in members
                j = colpos[f]
                @inbounds for t in (indptr[f] + 1):indptr[f + 1]
                    p = pos[ind[t - a + 1] + 1]
                    p == 0 || (out[p, j] = val[t - a + 1])
                end
            end
        end
    else
        # csr: compressed axis = cells (rows); indices are feature (column) ids
        fpos = position_lookup(ref.n_var, feature_idx)
        rowpos = position_lookup(ref.n_obs, cell_idx)
        runs = consecutive_runs(sort(cell_idx))
        _foreach_block(indptr, runs) do a, b, members
            ind = store_read(store, idx_path, a:b)
            val = store_read(store, data_path, a:b)
            for c in members
                i = rowpos[c]
                @inbounds for t in (indptr[c] + 1):indptr[c + 1]
                    p = fpos[ind[t - a + 1] + 1]
                    p == 0 || (out[i, p] = val[t - a + 1])
                end
            end
        end
    end
    return out
end

"""Sparse (cells × features) block from a CSR/CSC matrix element."""
function read_sparse_block_sparse(store::AbstractStore, ref::MatrixRef, indptr::Vector{Int},
                                  cell_idx::AbstractVector{Int}, feature_idx::AbstractVector{Int})
    I = Int[]; J = Int[]; V = Float64[]
    data_path = ref.path * "/data"; idx_path = ref.path * "/indices"
    if ref.kind == :csc
        pos = position_lookup(ref.n_obs, cell_idx)
        colpos = position_lookup(ref.n_var, feature_idx)
        _foreach_block(indptr, consecutive_runs(sort(feature_idx))) do a, b, members
            ind = store_read(store, idx_path, a:b); val = store_read(store, data_path, a:b)
            for f in members, t in (indptr[f] + 1):indptr[f + 1]
                p = pos[ind[t - a + 1] + 1]
                p == 0 || (push!(I, p); push!(J, colpos[f]); push!(V, val[t - a + 1]))
            end
        end
    else
        fpos = position_lookup(ref.n_var, feature_idx)
        rowpos = position_lookup(ref.n_obs, cell_idx)
        _foreach_block(indptr, consecutive_runs(sort(cell_idx))) do a, b, members
            ind = store_read(store, idx_path, a:b); val = store_read(store, data_path, a:b)
            for c in members, t in (indptr[c] + 1):indptr[c + 1]
                p = fpos[ind[t - a + 1] + 1]
                p == 0 || (push!(I, rowpos[c]); push!(J, p); push!(V, val[t - a + 1]))
            end
        end
    end
    return sparse(I, J, V, length(cell_idx), length(feature_idx))
end

"""Dense block from a dense matrix element stored python-shape (n_obs, n_var)."""
function read_dense_block(store::AbstractStore, ref::MatrixRef,
                          cell_idx::AbstractVector{Int}, feature_idx::AbstractVector{Int})
    out = Matrix{Float64}(undef, length(cell_idx), length(feature_idx))
    colpos = position_lookup(ref.n_var, feature_idx)
    for run in consecutive_runs(sort(feature_idx))
        # Julia dims are (n_var, n_obs): rows = features
        block = store_read(store, ref.path, run, :)           # length(run) × n_obs
        for (r, f) in enumerate(run)
            j = colpos[f]
            @inbounds for (i, c) in enumerate(cell_idx)
                out[i, j] = block[r, c]
            end
        end
    end
    return out
end

"""Dense block from a data-frame-encoded matrix (each feature is a column dataset)."""
function read_dataframe_block(store::AbstractStore, ref::MatrixRef,
                              cell_idx::AbstractVector{Int}, feature_idx::AbstractVector{Int})
    out = Matrix{Float64}(undef, length(cell_idx), length(feature_idx))
    for (j, f) in enumerate(feature_idx)
        col = read_element(store, ref.path * "/" * ref.columns[f])
        @inbounds for (i, c) in enumerate(cell_idx)
            out[i, j] = Float64(col[c])
        end
    end
    return out
end

"""Embedding block (cells × dims) from an obsm array stored python-shape (n_obs, k)."""
function read_array_block(store::AbstractStore, ref::ArrayRef, dims::AbstractVector{Int}, cell_idx::AbstractVector{Int})
    lo, hi = extrema(dims)
    block = store_read(store, ref.path, lo:hi, :)            # (hi-lo+1) × n_obs
    out = Matrix{Float64}(undef, length(cell_idx), length(dims))
    @inbounds for (j, d) in enumerate(dims), (i, c) in enumerate(cell_idx)
        out[i, j] = block[d - lo + 1, c]
    end
    return out
end
