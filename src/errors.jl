# Typed exceptions and the warning channel shared by every backend and wrapper.

abstract type SCUBAError <: Exception end

"""Thrown when none of the requested `vars` could be resolved."""
struct VarsNotFound <: SCUBAError
    missing::Vector{String}
    nrequested::Int
    suggestions::Dict{String,String}
end

"""Thrown by `fetch_metadata` when any requested metadata column is absent."""
struct MetadataNotFound <: SCUBAError
    missing::Vector{String}
    available::Vector{String}
end

struct LayerNotFound <: SCUBAError
    layer::String
    modality::String
    available::Vector{String}
end

struct ModalityNotFound <: SCUBAError
    modality::String
    available::Vector{String}
end

struct ReductionNotFound <: SCUBAError
    reduction::String
    available::Vector{String}
end

struct NoDefaultReduction <: SCUBAError
    available::Vector{String}
end

struct InvalidDims <: SCUBAError
    dims::Vector{Int}
    reduction::String
    ndims::Int
end

struct CellsNotFound <: SCUBAError
    missing::Vector{String}
    nrequested::Int
end

struct NonCategoricalVar <: SCUBAError
    var::String
    eltype::Type
end

struct ReadOnlySource <: SCUBAError
    source::String
end

struct NotImplementedBySource <: SCUBAError
    method::String
    source::String
end

_shown(v::AbstractVector{<:AbstractString}, n::Int=10) =
    length(v) > n ? " ($n out of $(length(v)) shown)" : ""
_join(v::AbstractVector{<:AbstractString}, n::Int=10) = join(first(v, n), ", ")

function Base.showerror(io::IO, e::VarsNotFound)
    head = length(e.missing) < e.nrequested ? "The following requested variables were not found" :
                                               "None of the requested variables were found"
    print(io, head, _shown(e.missing), ": ", _join(e.missing))
    if !isempty(e.suggestions)
        hints = ["'$k' -> '$v'" for (k, v) in e.suggestions]
        print(io, ". Did you mean: ", join(first(hints, 10), ", "), "?")
    end
end
Base.showerror(io::IO, e::MetadataNotFound) =
    print(io, "The following metadata variables were not found", _shown(e.missing), ": ",
          _join(e.missing), ". Use meta_varnames() to list available variables.")
Base.showerror(io::IO, e::LayerNotFound) =
    print(io, "Layer '", e.layer, "' does not exist in modality '", e.modality,
          "'. Available layers: ", join(e.available, ", "))
Base.showerror(io::IO, e::ModalityNotFound) =
    print(io, "Modality '", e.modality, "' not found. Available modalities: ", join(e.available, ", "))
Base.showerror(io::IO, e::ReductionNotFound) =
    print(io, "Reduction '", e.reduction, "' not found. Available reductions: ", join(e.available, ", "))
Base.showerror(io::IO, e::NoDefaultReduction) =
    print(io, "No default reduction (umap, tsne or pca) found in object. Available reductions: ",
          isempty(e.available) ? "none" : join(e.available, ", "))
Base.showerror(io::IO, e::InvalidDims) =
    print(io, "Dimensions ", e.dims, " are out of range for reduction '", e.reduction,
          "' which has ", e.ndims, " dimensions")
Base.showerror(io::IO, e::CellsNotFound) =
    print(io, "The following cells were not found in the object", _shown(e.missing), ": ", _join(e.missing))
Base.showerror(io::IO, e::NonCategoricalVar) =
    print(io, "Only categorical values permitted: metadata variable '", e.var,
          "' is not a categorical variable (eltype ", e.eltype, ")")
Base.showerror(io::IO, e::ReadOnlySource) = print(io, e.source, " does not support metadata assignment")
Base.showerror(io::IO, e::NotImplementedBySource) =
    print(io, e.method, " is not implemented for source type ", e.source)

# ---------------------------------------------------------------------------
# Warnings

"""
    SCUBAWarning(code, message, vars)

A structured warning. `code` is one of
`:duplicate_var`, `:some_vars_not_found`, `:alt_modality_renamed`,
`:ambiguous_var`, `:shadowed_feature`, `:layer_substituted`, `:large_result`.
"""
struct SCUBAWarning
    code::Symbol
    message::String
    vars::Vector{String}
end
SCUBAWarning(code::Symbol, message::AbstractString) = SCUBAWarning(code, String(message), String[])
SCUBAWarning(code::Symbol, message::AbstractString, var::AbstractString) =
    SCUBAWarning(code, String(message), [String(var)])

Base.show(io::IO, w::SCUBAWarning) = print(io, "SCUBAWarning(:", w.code, "): ", w.message)

const _WARNING_SINK = :scuba_warning_sink

"""
    scuba_warn(w::SCUBAWarning)

Emit a warning. Inside `with_collected_warnings` the warning is collected instead of
being logged, so host-language wrappers can re-emit it natively.
"""
function scuba_warn(w::SCUBAWarning)
    sink = get(task_local_storage(), _WARNING_SINK, nothing)
    if sink === nothing
        @warn w.message _group = :scuba code = w.code vars = w.vars
    else
        push!(sink::Vector{SCUBAWarning}, w)
    end
    return nothing
end
scuba_warn(code::Symbol, message::AbstractString, vars...) =
    scuba_warn(SCUBAWarning(code, String(message), String[string(v) for v in vars]))
scuba_warn(ws::AbstractVector{SCUBAWarning}) = (foreach(scuba_warn, ws); nothing)

"""
    with_collected_warnings(f) -> (result, Vector{SCUBAWarning})

Run `f()` and collect every SCUBA warning it emits instead of logging them.
"""
function with_collected_warnings(f)
    sink = SCUBAWarning[]
    result = task_local_storage(f, _WARNING_SINK, sink)
    return result, sink
end
