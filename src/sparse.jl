# Gather kernels shared by the in-memory, HDF5 and BPCells backends.

# Output sinks: a dense cells × features matrix, or one vector per feature.
@inline _put!(out::AbstractMatrix, i::Int, j::Int, v) = (@inbounds out[i, j] = v; nothing)
@inline _put!(cols::Vector{<:AbstractVector}, i::Int, j::Int, v) = (@inbounds cols[j][i] = v; nothing)
_columns(::Type{T}, ncells::Int, nfeat::Int) where {T} = Vector{T}[zeros(T, ncells) for _ in 1:nfeat]

"""
    SparseSink(T, nfeat)

Sparse output sink: one `(I, V)` buffer per feature column, filled through `_put!`. Explicit
zeros are dropped, so dense-backed layers yield genuinely sparse results and sparse and dense
results agree on values. Convert with [`to_sparse_matrix`](@ref) / [`to_sparse_columns`](@ref).
"""
struct SparseSink{T}
    I::Vector{Vector{Int}}
    V::Vector{Vector{T}}
end
SparseSink(::Type{T}, nfeat::Int) where {T} = SparseSink{T}([Int[] for _ in 1:nfeat], [T[] for _ in 1:nfeat])
@inline function _put!(s::SparseSink, i::Int, j::Int, v)
    iszero(v) && return nothing
    @inbounds (push!(s.I[j], i); push!(s.V[j], v))
    return nothing
end
# cell positions arrive in file order, i.e. already sorted unless `cells` was reordered
function _sorted_column(I::Vector{Int}, V::Vector{T}) where {T}
    issorted(I) && return I, V
    p = sortperm(I)
    return I[p], V[p]
end
"""`to_sparse_matrix(sink, ncells)` — cells × features `SparseMatrixCSC`."""
function to_sparse_matrix(s::SparseSink{T}, ncells::Int) where {T}
    nf = length(s.I)
    colptr = Vector{Int}(undef, nf + 1)
    colptr[1] = 1
    @inbounds for j in 1:nf
        colptr[j + 1] = colptr[j] + length(s.I[j])
    end
    nnz = colptr[end] - 1
    rowval = Vector{Int}(undef, nnz); nzval = Vector{T}(undef, nnz)
    @inbounds for j in 1:nf
        I, V = _sorted_column(s.I[j], s.V[j])
        isempty(I) && continue
        copyto!(rowval, colptr[j], I, 1, length(I))
        copyto!(nzval, colptr[j], V, 1, length(V))
    end
    return SparseMatrixCSC{T,Int}(ncells, nf, colptr, rowval, nzval)
end
"""`to_sparse_columns(sink, ncells)` — one `SparseVector` per feature."""
function _sparse_column(s::SparseSink{T}, j::Int, ncells::Int) where {T}
    I, V = _sorted_column(s.I[j], s.V[j])
    return SparseVector(ncells, I, V)
end
to_sparse_columns(s::SparseSink{T}, ncells::Int) where {T} = SparseVector{T,Int}[_sparse_column(s, j, ncells) for j in eachindex(s.I)]
"""Split a cells × features `SparseMatrixCSC` into its `SparseVector` columns (no copy of values per column beyond the slice)."""
_sparse_columns(X::SparseMatrixCSC{T}) where {T} = SparseVector{T,Int}[X[:, j] for j in 1:size(X, 2)]
"""Inverse of `_sparse_columns`: `SparseVector` columns -> `SparseMatrixCSC`."""
_hcat_sparse(cols::AbstractVector{<:SparseVector}) = isempty(cols) ? spzeros(0, 0) : reduce(hcat, cols)

"""
    RowMapSink(inner, rows)

Sink adaptor mapping the row (cell) position `i` of the writer to `rows[i]` in `inner`; used when
a modality holds only a subset of the cells (MuData).
"""
struct RowMapSink{S}
    inner::S
    rows::Vector{Int}
end
@inline _put!(s::RowMapSink, i::Int, j::Int, v) = _put!(s.inner, s.rows[i], j, v)

"""
    GroupSumSink(T, codes, ngroups, nfeat)

Reduction sink for [`aggregate`](@ref): `_put!(s, i, j, v)` adds `v` to `acc[codes[i], j]`, where
`codes[i]` is the 1-based group of the cell at position `i` of the request being written (`0`
skips the cell). The cells × features block is never materialised: every backend's `expression!`
streams straight into the `ngroups × nfeat` accumulator, which is all the memory a pseudobulk
profile needs. `GroupSumSink{T}(acc, codes)` wraps an existing accumulator so a request can be
written in several cell chunks (each with its own slice of `codes`).
"""
struct GroupSumSink{T}
    acc::Matrix{T}            # ngroups × nfeat
    codes::Vector{Int32}      # group per cell position of the current request
end
GroupSumSink(::Type{T}, codes::AbstractVector{<:Integer}, ngroups::Integer, nfeat::Integer) where {T} =
    GroupSumSink{T}(zeros(T, Int(ngroups), Int(nfeat)), Vector{Int32}(codes))
@inline function _put!(s::GroupSumSink, i::Int, j::Int, v)
    @inbounds g = s.codes[i]
    g == 0 && return nothing
    @inbounds s.acc[g, j] += v
    return nothing
end

"""`position_lookup(n, idx)` — vector `pos` with `pos[idx[k]] = k` and 0 elsewhere."""
function position_lookup(n::Integer, idx::AbstractVector{Int})
    pos = zeros(Int, n)
    @inbounds for (k, i) in enumerate(idx)
        pos[i] = k
    end
    return pos
end

"""Split a sorted vector of indices into runs of consecutive integers."""
function consecutive_runs(sorted_idx::AbstractVector{Int})
    runs = UnitRange{Int}[]
    isempty(sorted_idx) && return runs
    start = prev = sorted_idx[1]
    for i in Iterators.drop(sorted_idx, 1)
        if i == prev + 1
            prev = i
        else
            push!(runs, start:prev)
            start = prev = i
        end
    end
    push!(runs, start:prev)
    return runs
end

"""
    gather_cf(A, cell_idx, feature_idx) -> Matrix{Float64}

Dense block from a sparse matrix stored cells × features (AnnData orientation, columns are features).
"""
function gather_cf!(out, A::SparseMatrixCSC, cell_idx::AbstractVector{Int}, feature_idx::AbstractVector{Int})
    pos = position_lookup(size(A, 1), cell_idx)
    rows = rowvals(A); vals = nonzeros(A)
    @inbounds for (j, f) in enumerate(feature_idx)
        for k in nzrange(A, f)
            p = pos[rows[k]]
            p == 0 || _put!(out, p, j, vals[k])
        end
    end
    return out
end
gather_cf(A::SparseMatrixCSC, cell_idx::AbstractVector{Int}, feature_idx::AbstractVector{Int}) =
    gather_cf!(zeros(Float64, length(cell_idx), length(feature_idx)), A, cell_idx, feature_idx)

"""
    gather_fc(A, cell_idx, feature_idx) -> Matrix{Float64}

Dense block from a sparse matrix stored features × cells (Seurat/SCE/BPCells orientation, columns are cells).
"""
function gather_fc!(out, A::SparseMatrixCSC, cell_idx::AbstractVector{Int}, feature_idx::AbstractVector{Int})
    pos = position_lookup(size(A, 1), feature_idx)
    rows = rowvals(A); vals = nonzeros(A)
    @inbounds for (i, c) in enumerate(cell_idx)
        for k in nzrange(A, c)
            p = pos[rows[k]]
            p == 0 || _put!(out, i, p, vals[k])
        end
    end
    return out
end
gather_fc(A::SparseMatrixCSC, cell_idx::AbstractVector{Int}, feature_idx::AbstractVector{Int}) =
    gather_fc!(zeros(Float64, length(cell_idx), length(feature_idx)), A, cell_idx, feature_idx)

"""Sparse block (cells × features, values of type `T`) from a cells × features sparse matrix."""
gather_cf_sparse(A::SparseMatrixCSC, cell_idx::AbstractVector{Int}, feature_idx::AbstractVector{Int}, ::Type{T} = Float64) where {T} =
    to_sparse_matrix(gather_cf!(SparseSink(T, length(feature_idx)), A, cell_idx, feature_idx), length(cell_idx))

"""Sparse block (cells × features, values of type `T`) from a features × cells sparse matrix."""
gather_fc_sparse(A::SparseMatrixCSC, cell_idx::AbstractVector{Int}, feature_idx::AbstractVector{Int}, ::Type{T} = Float64) where {T} =
    to_sparse_matrix(gather_fc!(SparseSink(T, length(feature_idx)), A, cell_idx, feature_idx), length(cell_idx))

"""Dense block from any cells × features matrix."""
gather_dense(A::AbstractMatrix, cell_idx::AbstractVector{Int}, feature_idx::AbstractVector{Int}) =
    Matrix{Float64}(A[cell_idx, feature_idx])

"""Dense block from any features × cells matrix."""
gather_dense_t(A::AbstractMatrix, cell_idx::AbstractVector{Int}, feature_idx::AbstractVector{Int}) =
    Matrix{Float64}(permutedims(A[feature_idx, cell_idx]))

# Column (per-feature vector) variants
function gather_cells_by_features_columns(A::AbstractMatrix, cidx, fidx, ::Type{T}) where {T}
    cols = _columns(T, length(cidx), length(fidx))
    if A isa SparseMatrixCSC
        gather_cf!(cols, A, cidx, fidx)
    else
        @inbounds for (j, f) in enumerate(fidx), (i, c) in enumerate(cidx)
            cols[j][i] = A[c, f]
        end
    end
    return cols
end

# Convenience: dispatch on storage kind
gather_cells_by_features(A::SparseMatrixCSC, cidx, fidx) = gather_cf(A, cidx, fidx)
gather_cells_by_features(A::AbstractMatrix, cidx, fidx) = gather_dense(A, cidx, fidx)
gather_features_by_cells(A::SparseMatrixCSC, cidx, fidx) = gather_fc(A, cidx, fidx)
gather_features_by_cells(A::AbstractMatrix, cidx, fidx) = gather_dense_t(A, cidx, fidx)
