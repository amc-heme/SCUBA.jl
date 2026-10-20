# Native reader for BPCells matrix directories (no R or C++ dependency).
#
# A directory holds one sparse matrix with logical shape (nrow = features, ncol = cells),
# stored compressed along columns (`storage_order` = "col") or rows ("row"):
#   version        text: "(packed|unpacked)-(double|float|uint)-matrix-v(1|2)"
#   storage_order  text: "col" | "row"
#   shape          UINT32v1 header + [nrow, ncol]
#   row_names, col_names   newline-separated text (may be empty)
#   idxptr         UINT32v1 (v1) or UINT64v1 (v2) header + ncompressed+1 offsets
#   unpacked:      index (UINT32v1), val (DOUBLEv1 | FLOATSv1 | UINT32v1)
#   packed:        index_data, index_idx, index_starts (UINT32v1), index_idx_offsets (UINT64v1, v2 only)
#                  val: DOUBLEv1 | FLOATSv1 raw, or for uint: val_data, val_idx (UINT32v1), val_idx_offsets (UINT64v1)
# Every binary file starts with an 8-byte ASCII type header and is little-endian.

const BP_HEADERS = Dict(UInt32 => "UINT32v1", UInt64 => "UINT64v1", Float64 => "DOUBLEv1", Float32 => "FLOATSv1")

function _bp_mmap(path::AbstractString, ::Type{T}) where {T}
    isfile(path) || throw(ArgumentError("BPCells file not found: $path"))
    expected = BP_HEADERS[T]
    open(path, "r") do io
        hdr = String(read(io, 8))
        hdr == expected || error("BPCells file $path has header '$hdr', expected '$expected'")
        n = (filesize(path) - 8) ÷ sizeof(T)
        n == 0 && return T[]
        return Mmap.mmap(io, Vector{T}, n, 8)
    end
end

function _bp_text(path::AbstractString)
    isfile(path) || return String[]
    s = read(path, String)
    isempty(s) && return String[]
    lines = split(s, '\n'; keepempty = true)
    isempty(last(lines)) && pop!(lines)
    return String.(lines)
end

struct BP128Stream
    data::Vector{UInt32}          # packed words
    word_start::Vector{Int}       # 0-based word offset of each chunk (length nchunks+1)
    starts::Vector{UInt32}        # D1 chunk start values (empty for FOR)
    frame::UInt32                 # FOR frame (0 = plain)
    kind::Symbol                  # :d1 | :d1z | :for
end

function _bp128_stream(dir, prefix, kind; v2::Bool)
    data = _bp_mmap(joinpath(dir, prefix * "_data"), UInt32)
    idx = _bp_mmap(joinpath(dir, prefix * "_idx"), UInt32)
    offsets = v2 ? _bp_mmap(joinpath(dir, prefix * "_idx_offsets"), UInt64) : UInt64[0, typemax(UInt64)]
    ws = Vector{Int}(undef, length(idx))
    bump = 0
    boundaries = offsets[2:end]          # chunk numbers after which 2^32 must be added
    bi = 1
    for k in eachindex(idx)
        # entry k (1-based) is idx[k-1 (0-based)]; a boundary b applies to entries with 0-based index >= b
        while bi <= length(boundaries) && boundaries[bi] <= k - 1 && boundaries[bi] != typemax(UInt64)
            bump += 1; bi += 1
        end
        ws[k] = Int(idx[k]) + bump * (1 << 32)
    end
    starts = kind in (:d1, :d1z) ? _bp_mmap(joinpath(dir, prefix * "_starts"), UInt32) : UInt32[]
    return BP128Stream(data, ws, starts, kind == :for ? UInt32(1) : UInt32(0), kind)
end

nchunks(s::BP128Stream) = length(s.word_start) - 1

function decode_chunk!(out::Vector{UInt32}, s::BP128Stream, c::Int)   # c is 0-based chunk id
    woff = s.word_start[c + 1]
    bits = (s.word_start[c + 2] - woff) ÷ 4
    if s.kind == :d1
        bp128_unpack_d1!(out, s.data, woff, bits, s.starts[c + 1])
    elseif s.kind == :d1z
        bp128_unpack_d1z!(out, s.data, woff, bits, s.starts[c + 1])
    else
        bp128_unpack_for!(out, s.data, woff, bits, s.frame)
    end
    return out
end

"""
    BPCellsMatrix(dir) / open_bpcells(dir)

Lazy layer over a BPCells matrix directory. As a `LazyMatrix` it presents the cells × features
orientation expected by SCUBA (`size(m) == (ncells, nfeatures)`), while `row_names(m)` /
`col_names(m)` give the stored feature and cell names. Decoded index chunks are cached.
"""
mutable struct BPCellsMatrix <: LazyMatrix
    dir::String
    version::String
    packed::Bool
    transposed::Bool                 # storage_order == "row": compressed axis = logical rows
    nrow::Int                        # logical rows (features)
    ncol::Int                        # logical columns (cells)
    row_names::Vector{String}
    col_names::Vector{String}
    idxptr::Vector{Int}              # 0-based entry offsets, length ncompressed+1
    index::Union{Nothing,Vector{UInt32}}         # unpacked indices
    index_stream::Union{Nothing,BP128Stream}     # packed indices (d1)
    val::Union{Nothing,Vector{Float64},Vector{Float32},Vector{UInt32}}  # raw values
    val_stream::Union{Nothing,BP128Stream}       # packed uint values (FOR)
    chunk_cache::Dict{Int,Vector{UInt32}}
    val_chunk_cache::Dict{Int,Vector{UInt32}}
    lock::ReentrantLock
end

const BP_CHUNK_CACHE_MAX = Ref{Int}(8192)

function BPCellsMatrix(dir::AbstractString)
    isdir(dir) || throw(ArgumentError("BPCells directory not found: $dir"))
    version = strip(read(joinpath(dir, "version"), String))
    m = match(r"^(packed|unpacked)-(double|float|uint)-matrix-v(\d+)$", version)
    m === nothing && error("Unrecognised BPCells matrix version '$version' in $dir")
    packed = m[1] == "packed"
    vtype = m[2]
    v = parse(Int, m[3])
    order = strip(read(joinpath(dir, "storage_order"), String))
    order in ("col", "row") || error("storage_order must be 'col' or 'row', found '$order'")
    shape = _bp_mmap(joinpath(dir, "shape"), UInt32)
    nrow, ncol = Int(shape[1]), Int(shape[2])
    idxptr = v == 1 ? Vector{Int}(_bp_mmap(joinpath(dir, "idxptr"), UInt32)) :
                      Vector{Int}(_bp_mmap(joinpath(dir, "idxptr"), UInt64))
    rown = _bp_text(joinpath(dir, "row_names"))
    coln = _bp_text(joinpath(dir, "col_names"))
    index = nothing; istream = nothing; val = nothing; vstream = nothing
    if packed
        istream = _bp128_stream(dir, "index", :d1z; v2 = v >= 2)
        if vtype == "uint"
            vstream = _bp128_stream(dir, "val", :for; v2 = v >= 2)
        else
            val = _bp_mmap(joinpath(dir, "val"), vtype == "double" ? Float64 : Float32)
        end
    else
        index = _bp_mmap(joinpath(dir, "index"), UInt32)
        val = _bp_mmap(joinpath(dir, "val"), vtype == "double" ? Float64 : vtype == "float" ? Float32 : UInt32)
    end
    return BPCellsMatrix(String(dir), version, packed, order == "row", nrow, ncol, rown, coln, idxptr,
                         index, istream, val, vstream, Dict{Int,Vector{UInt32}}(), Dict{Int,Vector{UInt32}}(),
                         ReentrantLock())
end
open_bpcells(dir::AbstractString) = BPCellsMatrix(dir)

Base.size(m::BPCellsMatrix) = (m.ncol, m.nrow)            # cells × features
row_names(m::BPCellsMatrix) = m.row_names
col_names(m::BPCellsMatrix) = m.col_names
nnz_entries(m::BPCellsMatrix) = m.idxptr[end]
Base.show(io::IO, m::BPCellsMatrix) =
    print(io, "BPCellsMatrix(", m.nrow, " × ", m.ncol, ", ", m.version, ", ", m.transposed ? "row" : "col", "-major, ", nnz_entries(m), " nonzeros)")

function _cached_chunk(cache::Dict{Int,Vector{UInt32}}, stream::BP128Stream, c::Int, lk::ReentrantLock)
    lock(lk) do
        get!(cache, c) do
            length(cache) >= BP_CHUNK_CACHE_MAX[] && empty!(cache)
            decode_chunk!(Vector{UInt32}(undef, 128), stream, c)
        end
    end
end

# 0-based "other axis" index of entry t (0-based)
@inline function _entry_index(m::BPCellsMatrix, t::Int)
    m.index === nothing || return Int(m.index[t + 1])
    chunk = _cached_chunk(m.chunk_cache, m.index_stream, t >> 7, m.lock)
    return Int(chunk[(t & 127) + 1])
end
@inline function _entry_value(m::BPCellsMatrix, t::Int)
    m.val === nothing || return Float64(m.val[t + 1])
    chunk = _cached_chunk(m.val_chunk_cache, m.val_stream, t >> 7, m.lock)
    return Float64(chunk[(t & 127) + 1])
end

# Visit (other_index0, value) for every entry of compressed vector `comp` (1-based)
function _foreach_entry(f, m::BPCellsMatrix, comp::Int)
    a = m.idxptr[comp]; b = m.idxptr[comp + 1]
    for t in a:(b - 1)
        f(_entry_index(m, t), _entry_value(m, t))
    end
end

function pull(m::BPCellsMatrix, fidx::AbstractVector{Int}, cidx::AbstractVector{Int})
    out = zeros(Float64, length(cidx), length(fidx))
    if !m.transposed
        fpos = position_lookup(m.nrow, fidx)
        for (i, c) in enumerate(cidx)
            _foreach_entry(m, c) do f0, v
                p = fpos[f0 + 1]
                p == 0 || (out[i, p] = v)
            end
        end
    else
        cpos = position_lookup(m.ncol, cidx)
        for (j, f) in enumerate(fidx)
            _foreach_entry(m, f) do c0, v
                p = cpos[c0 + 1]
                p == 0 || (out[p, j] = v)
            end
        end
    end
    return out
end

function pull_sparse(m::BPCellsMatrix, fidx::AbstractVector{Int}, cidx::AbstractVector{Int})
    I = Int[]; J = Int[]; V = Float64[]
    if !m.transposed
        fpos = position_lookup(m.nrow, fidx)
        for (i, c) in enumerate(cidx)
            _foreach_entry(m, c) do f0, v
                p = fpos[f0 + 1]
                p == 0 || (push!(I, i); push!(J, p); push!(V, v))
            end
        end
    else
        cpos = position_lookup(m.ncol, cidx)
        for (j, f) in enumerate(fidx)
            _foreach_entry(m, f) do c0, v
                p = cpos[c0 + 1]
                p == 0 || (push!(I, p); push!(J, j); push!(V, v))
            end
        end
    end
    return sparse(I, J, V, length(cidx), length(fidx))
end

"""Read the whole matrix as a features × cells `SparseMatrixCSC` (BPCells orientation)."""
function SparseArrays.sparse(m::BPCellsMatrix)
    S = pull_sparse(m, collect(1:m.nrow), collect(1:m.ncol))
    return permutedims(S)
end

"""
    bpcells_modality(layers::AbstractDict{String,String}; features=nothing) -> (Modality, cells)

Build a `Modality` whose layers are BPCells directories (`layer name => dir`). Feature and
cell names are taken from the first directory and checked for agreement across layers.
"""
function bpcells_modality(layers::AbstractDict{<:AbstractString,<:AbstractString})
    mats = OrderedDict{String,LayerData}()
    feats = nothing; cells_ = nothing
    for (name, dir) in layers
        m = BPCellsMatrix(dir)
        if feats === nothing
            feats = copy(m.row_names); cells_ = copy(m.col_names)
        else
            m.row_names == feats || error("BPCells layer '$name' has different feature names")
            m.col_names == cells_ || error("BPCells layer '$name' has different cell names")
        end
        mats[String(name)] = m
    end
    return Modality(feats, mats), cells_
end
