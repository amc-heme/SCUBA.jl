# Live Python `anndata.AnnData` objects as SCUBA sources, through PythonCall.
# Slicing happens in Python (in memory or backed); only the requested block is converted.
module SCUBAPythonCallExt

using SCUBA, PythonCall, DataFrames, CategoricalArrays, SparseArrays, OrderedCollections
using SCUBA: AbstractSingleCellSource, SourceCache, ModalityNotFound, LayerNotFound, ReductionNotFound
import SCUBA: cell_ids, modality_keys, default_modality, default_layer, feature_names, layer_names,
              expression, metadata_names, metadata_column, reduction_keys, reduction_ndims, reduction,
              source_name, set_metadata!, invalidate!, pysource

const np = PythonCall.pynew()
const pd = PythonCall.pynew()
const sps = PythonCall.pynew()

function __init__()
    PythonCall.pycopy!(np, pyimport("numpy"))
    PythonCall.pycopy!(pd, pyimport("pandas"))
    PythonCall.pycopy!(sps, pyimport("scipy.sparse"))
    return nothing
end

struct PyModality
    key::String
    kind::Symbol                 # :X | :obsm | :raw
    features::Vector{String}
    layers::Vector{String}
end

"""
    PyAnnDataSource(adata::Py) / pysource(adata)

A live `anndata.AnnData` (in memory or backed) as a SCUBA source. Modalities: `"X"` (layers
`"X"` and `layers/*`), every `obsm` data frame, and `"raw"` when present. Reductions: every
numeric `obsm` array. Names are cached at construction; call `invalidate!` after mutating
the Python object's structure.
"""
mutable struct PyAnnDataSource <: AbstractSingleCellSource
    adata::Py
    cells::Vector{String}
    modalities::OrderedDict{String,PyModality}
    reductions::OrderedDict{String,Int}
    meta::Vector{String}
    cache::SourceCache
end

_strs(x) = String[pyconvert(String, pystr(e)) for e in pylist(x)]
_idx0(v::AbstractVector{Int}) = np.asarray(pylist(v .- 1), dtype = np.int64)

function PyAnnDataSource(adata::Py)
    pyhasattr(adata, "obs_names") && pyhasattr(adata, "var_names") ||
        throw(ArgumentError("object does not look like an anndata.AnnData"))
    cells = _strs(adata.obs_names.astype("str").to_list())
    var_names = _strs(adata.var_names.astype("str").to_list())
    mods = OrderedDict{String,PyModality}()
    xlayers = String["X"]
    append!(xlayers, _strs(adata.layers.keys()))
    mods["X"] = PyModality("X", :X, var_names, xlayers)
    reds = OrderedDict{String,Int}()
    for k in _strs(adata.obsm.keys())
        v = adata.obsm[k]
        if pyisinstance(v, pd.DataFrame)
            mods[k] = PyModality(k, :obsm, _strs(v.columns.astype("str").to_list()), ["X"])
        elseif pyhasattr(v, "shape") && pyconvert(Int, pylen(v.shape)) == 2
            reds[k] = pyconvert(Int, v.shape[1])
        end
    end
    if !pyis(adata.raw, pybuiltins.None)
        mods["raw"] = PyModality("raw", :raw, _strs(adata.raw.var_names.astype("str").to_list()), ["X"])
    end
    meta = _strs(adata.obs.columns.astype("str").to_list())
    return PyAnnDataSource(adata, cells, mods, reds, meta, SourceCache())
end
pysource(adata::Py) = PyAnnDataSource(adata)

source_name(::PyAnnDataSource) = "PyAnnDataSource"
cell_ids(s::PyAnnDataSource) = s.cells
modality_keys(s::PyAnnDataSource) = collect(keys(s.modalities))
default_modality(::PyAnnDataSource) = "X"
default_layer(::PyAnnDataSource) = "X"
function _mod(s::PyAnnDataSource, m)
    haskey(s.modalities, m) || throw(ModalityNotFound(String(m), modality_keys(s)))
    return s.modalities[m]
end
feature_names(s::PyAnnDataSource, m::AbstractString) = _mod(s, m).features
layer_names(s::PyAnnDataSource, m::AbstractString) = _mod(s, m).layers

function _to_dense(M::Py)
    if pyconvert(Bool, sps.issparse(M))
        M = M.toarray()
    elseif pyhasattr(M, "to_memory")           # backed dataset slices
        M = M.to_memory()
        pyconvert(Bool, sps.issparse(M)) && (M = M.toarray())
    end
    return pyconvert(Matrix{Float64}, np.ascontiguousarray(np.asarray(M, dtype = np.float64)))
end

function expression(s::PyAnnDataSource, m::AbstractString, layer::AbstractString,
                    fidx::AbstractVector{Int}, cidx::AbstractVector{Int})
    mod = _mod(s, m)
    layer in mod.layers || throw(LayerNotFound(String(layer), mod.key, mod.layers))
    ci = _idx0(collect(cidx)); fi = _idx0(collect(fidx))
    if mod.kind == :X
        sub = s.adata[ci, fi]
        M = layer == "X" ? sub.X : sub.layers[layer]
        return _to_dense(M)
    elseif mod.kind == :raw
        return _to_dense(s.adata.raw[ci, fi].X)
    else
        df = s.adata.obsm[mod.key]
        return _to_dense(df.iloc[ci, fi].to_numpy(dtype = np.float64))
    end
end

metadata_names(s::PyAnnDataSource) = s.meta

function _convert_series(col::Py)
    dtype = col.dtype
    if pyisinstance(dtype, pd.CategoricalDtype)
        codes = pyconvert(Vector{Int}, col.cat.codes.to_numpy())
        cats = _strs(col.cat.categories.astype("str").to_list())
        ordered = pyconvert(Bool, col.cat.ordered)
        return SCUBA._categorical_from_codes(codes, cats, ordered)
    end
    kind = pyconvert(String, pystr(dtype))
    if kind == "bool"
        return pyconvert(Vector{Bool}, col.to_numpy())
    elseif kind == "boolean"
        mask = pyconvert(Vector{Bool}, col.isna().to_numpy())
        vals = pyconvert(Vector{Bool}, col.fillna(false).to_numpy(dtype = np.bool_))
        return Union{Missing,Bool}[mask[i] ? missing : vals[i] for i in eachindex(vals)]
    elseif occursin(r"^(u?int|Int|UInt)", kind)
        if pyconvert(Bool, pd.api.types.is_extension_array_dtype(dtype))
            mask = pyconvert(Vector{Bool}, col.isna().to_numpy())
            vals = pyconvert(Vector{Int64}, col.fillna(0).to_numpy(dtype = np.int64))
            return Union{Missing,Int64}[mask[i] ? missing : vals[i] for i in eachindex(vals)]
        end
        return pyconvert(Vector{Int64}, col.to_numpy(dtype = np.int64))
    elseif startswith(kind, "float") || startswith(kind, "Float")
        if pyconvert(Bool, pd.api.types.is_extension_array_dtype(dtype))
            mask = pyconvert(Vector{Bool}, col.isna().to_numpy())
            vals = pyconvert(Vector{Float64}, col.fillna(0.0).to_numpy(dtype = np.float64))
            return Union{Missing,Float64}[mask[i] ? missing : vals[i] for i in eachindex(vals)]
        end
        return pyconvert(Vector{Float64}, col.to_numpy(dtype = np.float64))
    else
        mask = pyconvert(Vector{Bool}, col.isna().to_numpy())
        vals = _strs(col.astype("str").to_list())
        any(mask) || return vals
        return Union{Missing,String}[mask[i] ? missing : vals[i] for i in eachindex(vals)]
    end
end

function metadata_column(s::PyAnnDataSource, name::AbstractString, cidx::AbstractVector{Int})
    name in s.meta || throw(SCUBA.MetadataNotFound([String(name)], s.meta))
    col = s.adata.obs[name].iloc[_idx0(collect(cidx))]
    return _convert_series(col)
end

reduction_keys(s::PyAnnDataSource) = collect(keys(s.reductions))
function reduction_ndims(s::PyAnnDataSource, key::AbstractString)
    haskey(s.reductions, key) || throw(ReductionNotFound(String(key), reduction_keys(s)))
    return s.reductions[key]
end
function reduction(s::PyAnnDataSource, key::AbstractString, dims::AbstractVector{Int}, cidx::AbstractVector{Int})
    reduction_ndims(s, key)
    arr = np.asarray(s.adata.obsm[key], dtype = np.float64)
    block = arr[_idx0(collect(cidx))][pybuiltins.slice(nothing), _idx0(collect(dims))]
    return pyconvert(Matrix{Float64}, np.ascontiguousarray(block))
end

function set_metadata!(s::PyAnnDataSource, table)
    df = DataFrame(table; copycols = true)
    cc = DataFrames.metadata(df, SCUBA.CELL_COLUMN_KEY, nothing)
    cc === nothing || select!(df, Not(cc))
    nrow(df) == length(s.cells) || throw(DimensionMismatch("table has $(nrow(df)) rows; object has $(length(s.cells)) cells"))
    s.adata.obs = _to_pandas(df, s.cells)
    s.meta = _strs(s.adata.obs.columns.astype("str").to_list())
    invalidate!(s)
    return s
end

function invalidate!(s::PyAnnDataSource)
    empty!(s.cache)
    return s
end

# -----------------------------------------------------------------------------
# Julia DataFrame -> pandas (categoricals, missing and booleans preserved)

function _pyseries(v::AbstractVector)
    if v isa CategoricalVector
        codes = Int32[ismissing(x) ? Int32(-1) : Int32(levelcode(x) - 1) for x in v]
        cats = pylist(string.(levels(v)))
        return pd.Categorical.from_codes(np.asarray(pylist(codes), dtype = np.int32); categories = cats, ordered = isordered(v))
    end
    T = nonmissingtype(eltype(v))
    hasmissing = any(ismissing, v)
    if T <: Bool
        hasmissing || return np.asarray(pylist(Bool.(v)), dtype = np.bool_)
        return pd.array(pylist([ismissing(x) ? pybuiltins.None : Bool(x) for x in v]); dtype = "boolean")
    elseif T <: Integer
        hasmissing || return np.asarray(pylist(Int64.(v)), dtype = np.int64)
        return pd.array(pylist([ismissing(x) ? pybuiltins.None : Int64(x) for x in v]); dtype = "Int64")
    elseif T <: Real
        return np.asarray(pylist(Float64[ismissing(x) ? NaN : Float64(x) for x in v]), dtype = np.float64)
    else
        return pd.array(pylist([ismissing(x) ? pybuiltins.None : string(x) for x in v]); dtype = "object")
    end
end

function _to_pandas(df::DataFrame, cells::AbstractVector{<:AbstractString})
    d = pydict()
    for n in names(df)
        d[n] = _pyseries(df[!, n])
    end
    return pd.DataFrame(d; index = pd.Index(pylist(String.(cells))))
end

end # module
