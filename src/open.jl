# Open any supported container by sniffing its format.

"""
    open_source(path) -> source

Open an on-disk container with the right native reader: an `.h5ad` file (`H5ADSource`), a
MuData `.h5mu` file (`MuDataSource`), or a zarr directory (needs `using Zarr`). A BPCells matrix
directory is a single matrix, not a dataset: see `bpcells_modality`.
"""
function open_source(path::AbstractString)
    p = String(path)
    if isdir(p)
        if isfile(joinpath(p, "zarr.json")) || isfile(joinpath(p, ".zgroup"))
            hasmethod(open_zarr, Tuple{String}) || throw(ArgumentError("load Zarr.jl (`using Zarr`) to open zarr stores"))
            return open_zarr(p)
        end
        isfile(joinpath(p, "version")) && isfile(joinpath(p, "idxptr")) &&
            throw(ArgumentError("$p is a BPCells matrix directory (one matrix); build a dataset around it with bpcells_modality / open_bpcells"))
        throw(ArgumentError("unrecognised directory container: $p"))
    end
    isfile(p) || throw(ArgumentError("file not found: $p"))
    head = open(io -> read(io, 8), p)
    length(head) >= 6 && String(head[1:6]) == "MuData" && return open_h5mu(p)
    if head == UInt8[0x89, 0x48, 0x44, 0x46, 0x0d, 0x0a, 0x1a, 0x0a]
        file = h5open(p, "r")
        try
            store = H5Store(file)
            et = string(get(store_root_attrs(store), "encoding-type", ""))
            return et == "MuData" ? mudata_source_from_store(store, p) : anndata_source_from_store(store, p)
        catch
            close(file)
            rethrow()
        end
    end
    throw(ArgumentError("unrecognised file format: $p (expected .h5ad / .h5mu / zarr)"))
end
