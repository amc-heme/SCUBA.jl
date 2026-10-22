# Key table: how `{key}_{name}` prefixes map onto modalities, reductions and metadata.

struct KeyEntry
    kind::Symbol      # :modality | :reduction | :metadata
    target::String
end

const METADATA_KEY = "obs"

"""
    keytable(src) -> Dict{String,KeyEntry}

All addressable keys of a source: modality keys, reduction keys and the universal
metadata key `"obs"`. A modality key shadows a reduction key of the same name.
"""
function keytable(src::AbstractSingleCellSource)
    t = Dict{String,KeyEntry}()
    for m in modality_keys(src)
        t[String(m)] = KeyEntry(:modality, String(m))
    end
    for r in reduction_keys(src)
        haskey(t, r) || (t[String(r)] = KeyEntry(:reduction, String(r)))
    end
    for (a, target) in reduction_aliases(src)      # e.g. MuData: bare "X_umap" -> "rna_X_umap"
        haskey(t, a) || (t[String(a)] = KeyEntry(:reduction, String(target)))
    end
    haskey(t, METADATA_KEY) || (t[METADATA_KEY] = KeyEntry(:metadata, METADATA_KEY))
    return t
end

"""Keys whose `key_` prefix matches `var`, longest first."""
function matching_keys(kt::Dict{String,KeyEntry}, var::AbstractString)
    ks = String[k for k in keys(kt) if ncodeunits(var) > ncodeunits(k) + 1 &&
                 startswith(var, k) && codeunit(var, ncodeunits(k) + 1) == UInt8('_')]
    sort!(ks; by = ncodeunits, rev = true)
    return ks
end

"""Remainder of `var` after `key_`."""
strip_key(var::AbstractString, key::AbstractString) = String(SubString(var, ncodeunits(key) + 2))

"""
    all_keys(src) -> OrderedDict{String,String}

Modality and reduction names mapped to the key used to address them in `fetch_data`
(SCUBA.jl keys equal the names, so this is mostly for parity with the R API).
"""
function all_keys(src::AbstractSingleCellSource)
    d = OrderedDict{String,String}()
    for m in modality_keys(src)
        d[String(m)] = String(m)
    end
    for r in reduction_keys(src)
        haskey(d, r) || (d[String(r)] = String(r))
    end
    return d
end

# Case-insensitive suggestions for not-found vars, used in error messages.
function case_suggestions(src::AbstractSingleCellSource, missing::AbstractVector{<:AbstractString})
    isempty(missing) && return Dict{String,String}()
    lower = Dict{String,String}()
    for m in modality_keys(src), f in feature_names(src, m)
        get!(lower, lowercase(string(m, "_", f)), string(m, "_", f))
        m == default_modality(src) && get!(lower, lowercase(f), String(f))
    end
    for n in metadata_names(src)
        get!(lower, lowercase(n), String(n))
    end
    for r in reduction_keys(src), d in 1:reduction_ndims(src, r)
        get!(lower, lowercase(string(r, "_", d)), string(r, "_", d))
    end
    out = Dict{String,String}()
    for v in missing
        hit = get(lower, lowercase(v), nothing)
        hit === nothing || hit == v || (out[String(v)] = hit)
    end
    return out
end
