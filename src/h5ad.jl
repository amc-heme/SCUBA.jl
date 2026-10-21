# Backed .h5ad reader.

struct H5Modality
    features::Vector{String}
    layers::OrderedDict{String,MatrixRef}
end

"""
    H5ADSource(path) / open_h5ad(path)

Lazily read an `.h5ad` file. Only the schema (names, encodings, shapes) is read at open time;
expression values are read on demand, one coalesced hyperslab per block. Modalities: `"X"`
(layers `"X"` plus every `layers/*`), every `obsm` data frame (e.g. `protein`), and `raw` if
present. Reductions: every numeric `obsm` array (e.g. `X_umap`). Close with `close(src)`.
"""
mutable struct H5ADSource <: AbstractSingleCellSource
    path::String
    store::AbstractStore
    cells::Vector{String}
    obs_path::String
    obs_columns::Vector{String}
    modalities::OrderedDict{String,H5Modality}
    reductions::OrderedDict{String,ArrayRef}
    indptr_cache::Dict{String,Vector{Int}}
    cache::SourceCache
end

function H5ADSource(path::AbstractString)
    isfile(path) || throw(ArgumentError("file not found: $path"))
    file = h5open(String(path), "r")
    try
        return anndata_source_from_store(H5Store(file), String(path))
    catch
        close(file)
        rethrow()
    end
end

"""Build a source over any `AbstractStore` holding the anndata layout (used by the HDF5 and Zarr readers)."""
function anndata_source_from_store(store::AbstractStore, path::String)
    begin
        cells = dataframe_index(store, "obs")
        obs_cols = dataframe_columns(store, "obs")
        var_names = dataframe_index(store, "var")
        n_obs, n_var = length(cells), length(var_names)

        mods = OrderedDict{String,H5Modality}()
        xlayers = OrderedDict{String,MatrixRef}()
        store_haskey(store, "X") && (xlayers["X"] = matrix_ref(store, "X"))
        if store_haskey(store, "layers")
            for l in store_children(store, "layers")
                xlayers[l] = matrix_ref(store, "layers/" * l)
            end
        end
        mods["X"] = H5Modality(var_names, xlayers)

        reds = OrderedDict{String,ArrayRef}()
        if store_haskey(store, "obsm")
            for k in store_children(store, "obsm")
                p = "obsm/" * k
                et = encoding_type(store, p)
                if et == "dataframe"
                    ref = matrix_ref(store, p)
                    mods[k] = H5Modality(ref.columns, OrderedDict("X" => ref))
                elseif et == "array"
                    sz = store_size(store, p)
                    length(sz) == 2 || continue
                    reds[k] = ArrayRef(p, Int(sz[2]), Int(sz[1]))
                elseif et == "csr_matrix" || et == "csc_matrix"
                    ref = matrix_ref(store, p)   # sparse obsm: treat as modality with numbered features
                    mods[k] = H5Modality([string(i) for i in 1:ref.n_var], OrderedDict("X" => ref))
                end
            end
        end
        if store_haskey(store, "raw") && store_haskey(store, "raw/X")
            raw_var = dataframe_index(store, "raw/var")
            mods["raw"] = H5Modality(raw_var, OrderedDict("X" => matrix_ref(store, "raw/X")))
        end
        return H5ADSource(path, store, cells, "obs", obs_cols, mods, reds,
                          Dict{String,Vector{Int}}(), SourceCache())
    end
end
open_h5ad(path::AbstractString) = H5ADSource(path)

Base.close(src::H5ADSource) = (src.store isa H5Store && close(src.store.file); nothing)
Base.isopen(src::H5ADSource) = src.store isa H5Store ? isopen(src.store.file) : true

source_name(src::H5ADSource) = "H5ADSource($(basename(src.path)))"
cell_ids(src::H5ADSource) = src.cells
modality_keys(src::H5ADSource) = collect(keys(src.modalities))
default_modality(src::H5ADSource) = "X"
default_layer(src::H5ADSource) = "X"
function _h5mod(src::H5ADSource, m)
    haskey(src.modalities, m) || throw(ModalityNotFound(String(m), modality_keys(src)))
    return src.modalities[m]
end
feature_names(src::H5ADSource, m::AbstractString) = _h5mod(src, m).features
layer_names(src::H5ADSource, m::AbstractString) = collect(keys(_h5mod(src, m).layers))
function _h5ref(src::H5ADSource, m, layer)
    mod = _h5mod(src, m)
    haskey(mod.layers, layer) || throw(LayerNotFound(String(layer), String(m), collect(keys(mod.layers))))
    return mod.layers[layer]
end
function _indptr(src::H5ADSource, ref::MatrixRef)
    get!(src.indptr_cache, ref.path) do
        Vector{Int}(store_read(src.store, ref.path * "/indptr"))
    end
end
function expression(src::H5ADSource, m::AbstractString, layer::AbstractString,
                    fidx::AbstractVector{Int}, cidx::AbstractVector{Int})
    ref = _h5ref(src, m, layer)
    if ref.kind == :csr || ref.kind == :csc
        return read_sparse_block(src.store, ref, _indptr(src, ref), cidx, fidx)
    elseif ref.kind == :dense
        return read_dense_block(src.store, ref, cidx, fidx)
    else
        return read_dataframe_block(src.store, ref, cidx, fidx)
    end
end
function expression_columns(src::H5ADSource, m::AbstractString, layer::AbstractString,
                            fidx::AbstractVector{Int}, cidx::AbstractVector{Int}, ::Type{T}) where {T<:Real}
    ref = _h5ref(src, m, layer)
    cols = _columns(T, length(cidx), length(fidx))
    if ref.kind == :csr || ref.kind == :csc
        read_sparse_block!(cols, src.store, ref, _indptr(src, ref), cidx, fidx)
    elseif ref.kind == :dense
        read_dense_block!(cols, src.store, ref, cidx, fidx)
    else
        read_dataframe_block!(cols, src.store, ref, cidx, fidx)
    end
    return cols
end
supports_sparse(src::H5ADSource, m::AbstractString, layer::AbstractString) =
    _h5ref(src, m, layer).kind in (:csr, :csc)
function expression_sparse(src::H5ADSource, m::AbstractString, layer::AbstractString,
                           fidx::AbstractVector{Int}, cidx::AbstractVector{Int})
    ref = _h5ref(src, m, layer)
    ref.kind in (:csr, :csc) || return sparse(expression(src, m, layer, fidx, cidx))
    return read_sparse_block_sparse(src.store, ref, _indptr(src, ref), cidx, fidx)
end
metadata_names(src::H5ADSource) = src.obs_columns
function metadata_column(src::H5ADSource, name::AbstractString, cidx::AbstractVector{Int})
    name in src.obs_columns || throw(MetadataNotFound([String(name)], src.obs_columns))
    col = read_element(src.store, src.obs_path * "/" * name)
    return col[cidx]
end
reduction_keys(src::H5ADSource) = collect(keys(src.reductions))
function reduction_ndims(src::H5ADSource, key::AbstractString)
    haskey(src.reductions, key) || throw(ReductionNotFound(String(key), reduction_keys(src)))
    return src.reductions[key].k
end
function reduction(src::H5ADSource, key::AbstractString, dims::AbstractVector{Int}, cidx::AbstractVector{Int})
    haskey(src.reductions, key) || throw(ReductionNotFound(String(key), reduction_keys(src)))
    return read_array_block(src.store, src.reductions[key], dims, cidx)
end
