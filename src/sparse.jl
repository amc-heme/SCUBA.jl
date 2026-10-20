# Gather kernels shared by the in-memory, HDF5 and BPCells backends.

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
function gather_cf(A::SparseMatrixCSC, cell_idx::AbstractVector{Int}, feature_idx::AbstractVector{Int})
    out = zeros(Float64, length(cell_idx), length(feature_idx))
    pos = position_lookup(size(A, 1), cell_idx)
    rows = rowvals(A); vals = nonzeros(A)
    @inbounds for (j, f) in enumerate(feature_idx)
        for k in nzrange(A, f)
            p = pos[rows[k]]
            p == 0 || (out[p, j] = vals[k])
        end
    end
    return out
end

"""
    gather_fc(A, cell_idx, feature_idx) -> Matrix{Float64}

Dense block from a sparse matrix stored features × cells (Seurat/SCE/BPCells orientation, columns are cells).
"""
function gather_fc(A::SparseMatrixCSC, cell_idx::AbstractVector{Int}, feature_idx::AbstractVector{Int})
    out = zeros(Float64, length(cell_idx), length(feature_idx))
    pos = position_lookup(size(A, 1), feature_idx)
    rows = rowvals(A); vals = nonzeros(A)
    @inbounds for (i, c) in enumerate(cell_idx)
        for k in nzrange(A, c)
            p = pos[rows[k]]
            p == 0 || (out[i, p] = vals[k])
        end
    end
    return out
end

"""Sparse block (cells × features) from a cells × features sparse matrix."""
function gather_cf_sparse(A::SparseMatrixCSC, cell_idx::AbstractVector{Int}, feature_idx::AbstractVector{Int})
    pos = position_lookup(size(A, 1), cell_idx)
    rows = rowvals(A); vals = nonzeros(A)
    I = Int[]; J = Int[]; V = Float64[]
    @inbounds for (j, f) in enumerate(feature_idx)
        for k in nzrange(A, f)
            p = pos[rows[k]]
            p == 0 && continue
            push!(I, p); push!(J, j); push!(V, vals[k])
        end
    end
    return sparse(I, J, V, length(cell_idx), length(feature_idx))
end

"""Sparse block (cells × features) from a features × cells sparse matrix."""
function gather_fc_sparse(A::SparseMatrixCSC, cell_idx::AbstractVector{Int}, feature_idx::AbstractVector{Int})
    pos = position_lookup(size(A, 1), feature_idx)
    rows = rowvals(A); vals = nonzeros(A)
    I = Int[]; J = Int[]; V = Float64[]
    @inbounds for (i, c) in enumerate(cell_idx)
        for k in nzrange(A, c)
            p = pos[rows[k]]
            p == 0 && continue
            push!(I, i); push!(J, p); push!(V, vals[k])
        end
    end
    return sparse(I, J, V, length(cell_idx), length(feature_idx))
end

"""Dense block from any cells × features matrix."""
gather_dense(A::AbstractMatrix, cell_idx::AbstractVector{Int}, feature_idx::AbstractVector{Int}) =
    Matrix{Float64}(A[cell_idx, feature_idx])

"""Dense block from any features × cells matrix."""
gather_dense_t(A::AbstractMatrix, cell_idx::AbstractVector{Int}, feature_idx::AbstractVector{Int}) =
    Matrix{Float64}(permutedims(A[feature_idx, cell_idx]))

# Convenience: dispatch on storage kind
gather_cells_by_features(A::SparseMatrixCSC, cidx, fidx) = gather_cf(A, cidx, fidx)
gather_cells_by_features(A::AbstractMatrix, cidx, fidx) = gather_dense(A, cidx, fidx)
gather_features_by_cells(A::SparseMatrixCSC, cidx, fidx) = gather_fc(A, cidx, fidx)
gather_features_by_cells(A::AbstractMatrix, cidx, fidx) = gather_dense_t(A, cidx, fidx)
