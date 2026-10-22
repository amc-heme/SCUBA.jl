# AnnData zarr stores (zarr format 2 and 3) as SCUBA sources, through Zarr.jl.
# Only the store adapter lives here; all anndata decoding is shared with the HDF5 reader.
module SCUBAZarrExt

using SCUBA, Zarr, OrderedCollections
using SCUBA: AbstractStore, SourceCache, H5Modality, MatrixRef, ArrayRef, dataframe_index, dataframe_columns,
             matrix_ref, encoding_type, store_size
import SCUBA: store_haskey, store_isgroup, store_attrs, store_root_attrs, store_children, store_size, store_eltype, store_read, open_zarr

struct ZarrStore <: AbstractStore
    root::ZGroup
end

function _znav(s::ZarrStore, path::AbstractString)
    obj = s.root
    for part in split(path, '/'; keepempty = false)
        obj isa ZGroup || return nothing
        if haskey(obj.groups, part)
            obj = obj.groups[part]
        elseif haskey(obj.arrays, part)
            obj = obj.arrays[part]
        else
            return nothing
        end
    end
    return obj
end
function _zobj(s::ZarrStore, path::AbstractString)
    obj = _znav(s, path)
    obj === nothing && throw(KeyError(path))
    return obj
end

store_haskey(s::ZarrStore, path::AbstractString) = _znav(s, path) !== nothing
store_isgroup(s::ZarrStore, path::AbstractString) = _zobj(s, path) isa ZGroup
function store_children(s::ZarrStore, path::AbstractString)
    g = _zobj(s, path)
    g isa ZGroup || return String[]
    return vcat(collect(String, keys(g.groups)), collect(String, keys(g.arrays)))
end
store_size(s::ZarrStore, path::AbstractString) = size(_zobj(s, path))
store_eltype(s::ZarrStore, path::AbstractString) = eltype(_zobj(s, path))
function store_attrs(s::ZarrStore, path::AbstractString)
    a = _zobj(s, path).attrs
    out = Dict{String,Any}()
    for (k, v) in a
        out[String(k)] = v isa AbstractVector ? collect(v) : v
    end
    return out
end
function store_read(s::ZarrStore, path::AbstractString)
    arr = _zobj(s, path)
    return arr[ntuple(_ -> Colon(), ndims(arr))...]
end
store_read(s::ZarrStore, path::AbstractString, idx...) = _zobj(s, path)[idx...]

"""
    open_zarr(path) -> H5ADSource-like source over a zarr store

Open an AnnData zarr store (format 2 or 3). Returns a `SCUBA.H5ADSource` whose store is a
`ZarrStore`, so all accessors behave exactly as for `.h5ad` files.
"""
function open_zarr(path::AbstractString)
    root = zopen(String(path), "r")
    root isa ZGroup || throw(ArgumentError("$path is not a zarr group"))
    store = ZarrStore(root)
    et = string(get(root.attrs, "encoding-type", ""))
    return et == "MuData" ? SCUBA.mudata_source_from_store(store, String(path)) :
                            SCUBA.anndata_source_from_store(store, String(path))
end
store_root_attrs(s::ZarrStore) = Dict{String,Any}(s.root.attrs)

end # module
