# Chunked iteration over cells for results too wide to hold at once.

"""Default byte budget per chunk for [`eachchunk`](@ref) (256 MiB)."""
const CHUNK_BYTES = Ref{Int}(256 * 2^20)

struct ChunkIterator{S<:AbstractSingleCellSource}
    src::S
    found::Vector{ResolvedVar}
    cell_idx::Vector{Int}
    cell_names::Vector{String}
    fixed::Dict{Int,AbstractVector}     # metadata / reduction columns over all requested cells, by position
    n::Int                              # cells per chunk
    eltype::DataType
    sparse::Bool
end

"""
    eachchunk(src, vars; layer=nothing, cells=nothing, strict=false, eltype=Float64, sparse=false,
              cells_per_chunk=nothing, chunk_bytes=CHUNK_BYTES[]) -> iterator of DataFrames

Iterate the result of `fetch_data(src, vars; ...)` in consecutive slices of the requested cell
order, so a result far wider than memory can be consumed one block at a time. Every chunk has
the same columns, types and metadata as the full table, so `reduce(vcat, collect(it))` equals
`fetch_data(src, vars; ...)` (without the dense-size warning). Variables are resolved once and
resolution warnings are emitted once, by this call; metadata and reduction columns are read once
and sliced. The chunk size is `cells_per_chunk` or derived from `chunk_bytes` and the number of
expression columns.
"""
function eachchunk(src::AbstractSingleCellSource, vars::AbstractVector{<:AbstractString};
                   layer = nothing, cells = nothing, strict::Bool = false, eltype::Type{<:Real} = Float64,
                   sparse::Bool = false, cells_per_chunk::Union{Nothing,Integer} = nothing,
                   chunk_bytes::Integer = CHUNK_BYTES[])
    _check_eltype(src, eltype, sparse)
    res = resolve_vars(src, vars; layer)
    if strict && !isempty(res.missing)
        throw(VarsNotFound(res.missing, length(vars), case_suggestions(src, res.missing)))
    end
    cell_idx, cell_names = resolve_cells(src, cells)
    scuba_warn(res.warnings)
    found = res.found
    nexpr = count(rv -> rv.kind == :expression, found)
    nred = count(rv -> rv.kind == :reduction, found)
    bytes_per_cell = max(1, sizeof(eltype) * nexpr + 8 * nred + 8)
    n = cells_per_chunk === nothing ? clamp(Int(chunk_bytes) ÷ bytes_per_cell, 1, max(1, length(cell_idx))) :
                                      max(1, Int(cells_per_chunk))
    fixed = Dict{Int,AbstractVector}()
    red_groups = OrderedDict{String,Vector{Int}}()
    for (p, rv) in enumerate(found)
        if rv.kind == :metadata
            fixed[p] = metadata_column(src, rv.target, cell_idx)
        elseif rv.kind == :reduction
            push!(get!(red_groups, rv.target, Int[]), p)
        end
    end
    for (key, ps) in red_groups
        M = reduction(src, key, [found[p].index for p in ps], cell_idx)
        for (j, p) in enumerate(ps)
            fixed[p] = M[:, j]
        end
    end
    return ChunkIterator(src, collect(found), collect(cell_idx), collect(String, cell_names), fixed, n, eltype, sparse)
end
eachchunk(src::AbstractSingleCellSource, var::AbstractString; kw...) = eachchunk(src, [var]; kw...)

Base.length(it::ChunkIterator) = max(1, cld(length(it.cell_idx), it.n))
Base.eltype(::Type{<:ChunkIterator}) = DataFrame
function Base.iterate(it::ChunkIterator, start::Int = 1)
    N = length(it.cell_idx)
    start > max(N, 1) && return nothing
    stop = min(start + it.n - 1, N)
    range = start:stop
    pre = Dict{Int,AbstractVector}(p => v[range] for (p, v) in it.fixed)
    df = _fetch_resolved(it.src, it.found, it.cell_idx[range], it.cell_names[range], it.eltype, it.sparse; prefetched = pre)
    return df, stop + 1
end
