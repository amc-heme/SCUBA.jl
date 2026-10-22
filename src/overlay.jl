# A metadata overlay over any source: the inner source serves expression and reductions, the
# overlay owns the metadata table. Gives mutable metadata (`set_metadata!`) to read-only file
# sources such as `H5ADSource`, and lets a host (the R package) keep the metadata on its own side
# by handing the overlay an empty table so metadata names are never resolved in Julia.

"""
    MetadataOverlay(inner; obs = metadata_table(inner, 1:ncells(inner)))

Source that forwards everything to `inner` except metadata, which comes from the `obs`
`DataFrame` (one row per cell of `inner`). `set_metadata!` replaces the table.
"""
mutable struct MetadataOverlay{S<:AbstractSingleCellSource} <: AbstractSingleCellSource
    inner::S
    obs::DataFrame
    function MetadataOverlay(inner::S, obs::DataFrame) where {S<:AbstractSingleCellSource}
        n = ncells(inner)
        (ncol(obs) == 0 || nrow(obs) == n) ||
            throw(DimensionMismatch("overlay table has $(nrow(obs)) rows; source has $n cells"))
        return new{S}(inner, obs)
    end
end
MetadataOverlay(inner::AbstractSingleCellSource; obs::Union{Nothing,DataFrame} = nothing) =
    MetadataOverlay(inner, obs === nothing ? metadata_table(inner, 1:ncells(inner)) : DataFrame(obs; copycols = true))

# forwarded
cell_ids(o::MetadataOverlay) = cell_ids(o.inner)
modality_keys(o::MetadataOverlay) = modality_keys(o.inner)
default_modality(o::MetadataOverlay) = default_modality(o.inner)
default_layer(o::MetadataOverlay) = default_layer(o.inner)
feature_names(o::MetadataOverlay, m::AbstractString) = feature_names(o.inner, m)
layer_names(o::MetadataOverlay, m::AbstractString) = layer_names(o.inner, m)
layer_exists(o::MetadataOverlay, m::AbstractString, l::AbstractString) = layer_exists(o.inner, m, l)
expression(o::MetadataOverlay, m::AbstractString, l::AbstractString, fidx::AbstractVector{Int}, cidx::AbstractVector{Int}) =
    expression(o.inner, m, l, fidx, cidx)
expression_columns(o::MetadataOverlay, m::AbstractString, l::AbstractString, fidx::AbstractVector{Int},
                   cidx::AbstractVector{Int}, ::Type{T}) where {T<:Real} = expression_columns(o.inner, m, l, fidx, cidx, T)
supports_sparse(o::MetadataOverlay, m::AbstractString, l::AbstractString) = supports_sparse(o.inner, m, l)
expression_sparse(o::MetadataOverlay, m::AbstractString, l::AbstractString, fidx::AbstractVector{Int}, cidx::AbstractVector{Int}) =
    expression_sparse(o.inner, m, l, fidx, cidx)
expression_sparse(o::MetadataOverlay, m::AbstractString, l::AbstractString, fidx::AbstractVector{Int}, cidx::AbstractVector{Int}, ::Type{T}) where {T} =
    expression_sparse(o.inner, m, l, fidx, cidx, T)
expression!(out, o::MetadataOverlay, m::AbstractString, l::AbstractString, fidx::AbstractVector{Int}, cidx::AbstractVector{Int}) =
    expression!(out, o.inner, m, l, fidx, cidx)
reduction_keys(o::MetadataOverlay) = reduction_keys(o.inner)
reduction_ndims(o::MetadataOverlay, k::AbstractString) = reduction_ndims(o.inner, k)
reduction(o::MetadataOverlay, k::AbstractString, dims::AbstractVector{Int}, cidx::AbstractVector{Int}) = reduction(o.inner, k, dims, cidx)
reduction_dimnames_impl(o::MetadataOverlay, k::AbstractString, dims::AbstractVector{Int}) = reduction_dimnames_impl(o.inner, k, dims)
reduction_search_names(o::MetadataOverlay, k::AbstractString) = reduction_search_names(o.inner, k)
reduction_aliases(o::MetadataOverlay) = reduction_aliases(o.inner)
check_eltype(o::MetadataOverlay, ::Type{T}) where {T} = check_eltype(o.inner, T)
canonical_modality(o::MetadataOverlay, name::AbstractString) = canonical_modality(o.inner, name)
canonical_reduction(o::MetadataOverlay, name::AbstractString) = canonical_reduction(o.inner, name)
cache(o::MetadataOverlay) = cache(o.inner)
source_name(o::MetadataOverlay) = "MetadataOverlay(" * source_name(o.inner) * ")"
Base.close(o::MetadataOverlay) = (hasmethod(close, Tuple{typeof(o.inner)}) && close(o.inner); nothing)
Base.isopen(o::MetadataOverlay) = hasmethod(isopen, Tuple{typeof(o.inner)}) ? isopen(o.inner) : true

# owned
metadata_names(o::MetadataOverlay) = names(o.obs)
function metadata_column(o::MetadataOverlay, name::AbstractString, cidx::AbstractVector{Int})
    hasproperty(o.obs, name) || throw(MetadataNotFound([String(name)], names(o.obs)))
    return o.obs[cidx, name]
end
metadata_table(o::MetadataOverlay, cidx::AbstractVector{Int}) = o.obs[cidx, :]
function set_metadata!(o::MetadataOverlay, table)
    df = DataFrame(table; copycols = true)
    cc = metadata(df, CELL_COLUMN_KEY, nothing)
    cc === nothing || select!(df, Not(cc))
    nrow(df) == ncells(o) || throw(DimensionMismatch("table has $(nrow(df)) rows; object has $(ncells(o)) cells"))
    o.obs = df
    invalidate!(o)
    return o
end
invalidate!(o::MetadataOverlay) = (invalidate!(o.inner); o)
