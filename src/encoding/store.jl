# Minimal hierarchical-store API used by the AnnData decoders, so the same decoding code
# serves HDF5 (.h5ad) and, through an extension, Zarr stores.

abstract type AbstractStore end

"""`store_haskey(store, path)`"""
function store_haskey end
"""`store_isgroup(store, path)`"""
function store_isgroup end
"""`store_attrs(store, path) -> Dict{String,Any}`"""
function store_attrs end
"""`store_children(store, path) -> Vector{String}`"""
function store_children end
"""`store_size(store, path) -> Tuple` (Julia dimension order)"""
function store_size end
"""`store_read(store, path)` — whole dataset"""
function store_read end
"""`store_read(store, path, idx...)` — slice with ranges/colons in Julia dimension order"""
function store_eltype end

# ---------------------------------------------------------------------------
# HDF5 implementation

struct H5Store <: AbstractStore
    file::HDF5.File
end

_h5obj(s::H5Store, path::AbstractString) = s.file[String(path)]
store_haskey(s::H5Store, path::AbstractString) = haskey(s.file, String(path))
store_isgroup(s::H5Store, path::AbstractString) = _h5obj(s, path) isa HDF5.Group
store_children(s::H5Store, path::AbstractString) = collect(String, keys(_h5obj(s, path)))
store_size(s::H5Store, path::AbstractString) = size(_h5obj(s, path))
store_eltype(s::H5Store, path::AbstractString) = eltype(_h5obj(s, path))

function store_attrs(s::H5Store, path::AbstractString)
    obj = _h5obj(s, path)
    out = Dict{String,Any}()
    for k in keys(HDF5.attrs(obj))
        out[String(k)] = _normalize_attr(HDF5.read_attribute(obj, k))
    end
    return out
end
_normalize_attr(x) = x
_normalize_attr(x::AbstractArray{<:AbstractString}) = collect(String, x)
_normalize_attr(x::AbstractString) = String(x)

# h5py stores booleans as a 2-member enum over int8; HDF5.jl reads them as Int8.
function _is_h5py_bool(ds)
    ds isa HDF5.Dataset || return false
    dt = HDF5.datatype(ds)
    try
        HDF5.API.h5t_get_class(dt) == HDF5.API.H5T_ENUM || return false
        HDF5.API.h5t_get_nmembers(dt) == 2 || return false
        names = Set(uppercase(HDF5.API.h5t_get_member_name(dt, i)) for i in 0:1)
        return names == Set(["FALSE", "TRUE"])
    finally
        close(dt)
    end
end
_boolify(ds, x) = _is_h5py_bool(ds) ? (x isa AbstractArray ? Bool.(x) : Bool(x)) : x

function store_read(s::H5Store, path::AbstractString)
    ds = _h5obj(s, path)
    return _normalize_data(_boolify(ds, read(ds)))
end
function store_read(s::H5Store, path::AbstractString, idx...)
    ds = _h5obj(s, path)
    return _normalize_data(_boolify(ds, ds[idx...]))
end

_normalize_data(x) = x
_normalize_data(x::AbstractArray{<:AbstractString}) = collect(String, x)
_normalize_data(x::AbstractString) = String(x)
