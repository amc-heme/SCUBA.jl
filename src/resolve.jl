# Variable resolution: turns user-facing `vars` into typed fetch instructions.
# This is the single implementation of the rules that the R package spreads over
# SeuratObject::FetchData, fetch_data.SingleCellExperiment and fetch_anndata.py.

struct ResolvedVar
    input::String      # var as written by the caller
    output::String     # column name in the result
    kind::Symbol       # :expression | :reduction | :metadata
    target::String     # modality key, reduction key, or metadata column name
    layer::Union{Nothing,String}
    index::Int         # feature index or reduction dim (0 for metadata)
end

struct ResolveResult
    found::Vector{ResolvedVar}
    missing::Vector{String}
    warnings::Vector{SCUBAWarning}
end

_dedup_key(rv::ResolvedVar) = (rv.kind, rv.target, rv.layer, rv.index)

"""
    resolve_vars(src, vars; layer=nothing) -> ResolveResult

Classify each var as a keyed or unkeyed feature, a reduction dimension or a metadata
column, following the unified SCUBA rules:

1. A key matches only as `key_` (with the underscore); matching is exact and case-sensitive.
2. Nested keys (`X` vs `X_umap`) are tried longest-first; the first whose remainder exists wins.
3. Unkeyed names resolve metadata → main modality → a single alternate modality (renamed
   `key_var` with a warning); names present in several alternates are ambiguous.
4. Duplicates are collapsed by resolved target, not by spelling.

Throws `VarsNotFound` when nothing resolves. Warnings are returned in the result (not
emitted); `fetch_data` emits them.
"""
function resolve_vars(src::AbstractSingleCellSource, vars::AbstractVector{<:AbstractString}; layer = nothing)
    W = SCUBAWarning[]
    kt = keytable(src)
    main = default_modality(src)
    meta = metadata_names(src)
    metaset = Set{String}(String.(meta))
    found = ResolvedVar[]
    missing_vars = String[]
    seen = Set{Tuple{Symbol,String,Union{Nothing,String},Int}}()
    layer_cache = Dict{String,String}()
    layer_for(mod) = get!(layer_cache, mod) do
        resolve_layer_collect!(W, src, mod, layer)
    end

    # literal duplicates
    literal_dups = String[]
    uniq = String[]
    seen_literal = Set{String}()
    for v in vars
        s = String(v)
        s in seen_literal ? push!(literal_dups, s) : (push!(seen_literal, s); push!(uniq, s))
    end
    isempty(literal_dups) || push!(W, SCUBAWarning(:duplicate_var,
        "The following entries to `vars` are duplicates: '" * join(unique(literal_dups), "', '") *
        "'. Only one entry for each duplicated var will be returned.", unique(literal_dups)))

    function push_found!(rv::ResolvedVar)
        k = _dedup_key(rv)
        if k in seen
            push!(W, SCUBAWarning(:duplicate_var,
                "'$(rv.input)' refers to the same data as an earlier entry and was dropped.", rv.input))
        else
            push!(seen, k)
            push!(found, rv)
        end
        return nothing
    end

    for v in uniq
        # (a) keyed
        hits = ResolvedVar[]
        cands = matching_keys(kt, v)
        for k in cands
            rest = strip_key(v, k)
            e = kt[k]
            if e.kind == :modality
                i = get(feature_index(src, e.target), rest, 0)
                i == 0 && continue
                push!(hits, ResolvedVar(v, v, :expression, e.target, layer_for(e.target), i))
            elseif e.kind == :reduction
                d = tryparse(Int, rest)
                (d === nothing || d < 1 || d > reduction_ndims(src, e.target)) && continue
                push!(hits, ResolvedVar(v, only(reduction_dimnames_impl(src, e.target, [d])), :reduction, e.target, nothing, d))
            else
                rest in metaset || continue
                push!(hits, ResolvedVar(v, rest, :metadata, rest, nothing, 0))
            end
        end
        if !isempty(hits)
            if length(hits) > 1
                others = join(["$(h.kind) via key '$(cands[j])'" for (j, h) in enumerate(hits) if j > 1], ", ")
                push!(W, SCUBAWarning(:ambiguous_var,
                    "'$v' matches several keys; resolved as $(hits[1].kind) via key '$(cands[1])' (also: $others).", v))
            end
            push_found!(hits[1])
            continue
        end
        # (b) unkeyed metadata
        if v in metaset
            if haskey(feature_index(src, main), v)
                push!(W, SCUBAWarning(:shadowed_feature,
                    "'$v' is both a metadata variable and a feature in modality '$main'. The metadata " *
                    "variable was returned; use '$(main)_$v' to fetch the feature.", v))
            end
            push_found!(ResolvedVar(v, v, :metadata, v, nothing, 0))
            continue
        end
        # (c) unkeyed main modality feature
        i = get(feature_index(src, main), v, 0)
        if i != 0
            push_found!(ResolvedVar(v, v, :expression, main, layer_for(main), i))
            continue
        end
        # (d) unkeyed alternate modality feature
        alts = Tuple{String,Int}[]
        for m in modality_keys(src)
            m == main && continue
            j = get(feature_index(src, m), v, 0)
            j == 0 || push!(alts, (String(m), j))
        end
        if length(alts) == 1
            m, j = alts[1]
            out = string(m, "_", v)
            push!(W, SCUBAWarning(:alt_modality_renamed,
                "$v was passed to fetch_data without specifying the key of the modality to pull the " *
                "feature from, and it is not present in the main modality. The feature was found in " *
                "'$m' and successfully returned as '$out'.", v))
            push_found!(ResolvedVar(v, out, :expression, m, layer_for(m), j))
            continue
        elseif length(alts) > 1
            push!(W, SCUBAWarning(:ambiguous_var,
                "'$v' was found in multiple modalities (" * join(first.(alts), ", ") *
                ") and will not be returned. Specify a key, e.g. '" * string(alts[1][1], "_", v) * "'.", v))
        end
        push!(missing_vars, v)
    end

    if isempty(found)
        throw(VarsNotFound(missing_vars, length(vars), case_suggestions(src, missing_vars)))
    end
    if !isempty(missing_vars)
        push!(W, SCUBAWarning(:some_vars_not_found,
            "The following requested variables were not found" * _shown(missing_vars) * ": " *
            _join(missing_vars), missing_vars))
    end
    return ResolveResult(found, missing_vars, W)
end

resolve_vars(src::AbstractSingleCellSource, var::AbstractString; kw...) = resolve_vars(src, [var]; kw...)

# resolve_layer, but collecting the substitution warning into W instead of emitting it
function resolve_layer_collect!(W::Vector{SCUBAWarning}, src, modality, layer)
    result, ws = with_collected_warnings() do
        resolve_layer(src, modality, layer)
    end
    append!(W, ws)
    return result
end
