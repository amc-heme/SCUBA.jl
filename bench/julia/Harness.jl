# Shared helpers for the Julia benchmark harness: configuration, dataset metadata, timing,
# memory, checksums and the result-row CSV writer (schema in common/schema.md).
module Harness

using SCUBA, DataFrames, CategoricalArrays, SparseArrays, SHA, JSON, TOML, Dates, Statistics

export Ctx, load_ctx, result_row!, timed, cells_for_subset, feature_checksum, metadata_checksum,
       reduction_checksum, fetch_cells_checksum, unique_values_checksum, checksum_ok, skip_rule_matches,
       existing_keys, row_key, cpu_time_ms, rss_mb, maxrss_mb, bench_dir, data_dir

const SCHEMA = ["suite_version", "run_id", "timestamp", "host", "language", "implementation", "engine", "storage",
                "scale", "workload", "n_features", "gene_set", "cell_subset", "layer", "n_cells_out", "cache_state",
                "phase", "rep", "threads", "wall_ms", "cpu_ms", "open_ms", "rss_mb", "maxrss_mb", "alloc_mb", "gc_ms",
                "result_checksum", "checksum_ok", "status", "note"]

bench_dir() = normpath(joinpath(@__DIR__, ".."))

mutable struct Ctx
    cfg::Dict{String,Any}
    scale::String
    run_id::String
    data::String                       # dataset directory
    selection::Dict{String,Any}
    expected::Dict{String,Any}
    manifest::Dict{String,Any}
    cache_state::String
    out_dir::String
    host::String
    io::Union{Nothing,IO}
    out_path::String
    existing::Set{String}
end

data_dir(cfg) = expanduser(cfg["suite"]["data_dir"])

function load_ctx(config_path, scale; run_id = "", cache_state = "warm", out_name = "")
    cfg = TOML.parsefile(config_path)
    d = joinpath(data_dir(cfg), scale)
    man = JSON.parsefile(joinpath(d, "MANIFEST.json"))
    man["complete"] == true || error("dataset $scale is incomplete")
    sel = JSON.parsefile(joinpath(d, "selection.json"))
    exp = JSON.parsefile(joinpath(d, "expected.json"))
    rid = isempty(run_id) ? (isempty(cfg["suite"]["run_id"]) ? Dates.format(now(UTC), "yyyymmddTHHMMSS") : cfg["suite"]["run_id"]) : run_id
    out_dir = joinpath(bench_dir(), cfg["suite"]["results_dir"], "raw", rid)
    mkpath(out_dir)
    ctx = Ctx(cfg, scale, rid, d, sel, exp, man, cache_state, out_dir, gethostname(), nothing, "", Set{String}())
    isempty(out_name) || open_output!(ctx, out_name)
    return ctx
end

function open_output!(ctx::Ctx, name::AbstractString)
    ctx.out_path = joinpath(ctx.out_dir, name * ".csv")
    fresh = !isfile(ctx.out_path) || filesize(ctx.out_path) == 0
    ctx.existing = fresh ? Set{String}() : existing_keys(ctx.out_path)
    ctx.io = open(ctx.out_path, "a")
    fresh && println(ctx.io, join(SCHEMA, ","))
    flush(ctx.io)
    return ctx
end

Base.close(ctx::Ctx) = (ctx.io === nothing || close(ctx.io); ctx.io = nothing)

_csv(x) = (s = x === nothing || x === missing ? "" : string(x); occursin(r"[,\"\n]", s) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s)
_fmt(x::AbstractFloat) = isnan(x) ? "" : string(round(x; digits = 4))
_fmt(x) = x

"""Key identifying a measurement cell (used for --resume)."""
row_key(implementation, workload, n_features, gene_set, cell_subset, layer, cache_state, phase, rep) =
    join((implementation, workload, n_features, gene_set, cell_subset, layer, cache_state, phase, rep), "|")

function existing_keys(path)
    keys = Set{String}()
    isfile(path) || return keys
    hdr = nothing
    for (i, line) in enumerate(eachline(path))
        if i == 1
            hdr = split(line, ",")
            continue
        end
        f = split(line, ",")              # keys never contain commas
        length(f) >= length(hdr) || continue
        col(n) = f[findfirst(==(n), hdr)]
        push!(keys, row_key(col("implementation"), col("workload"), col("n_features"), col("gene_set"),
                            col("cell_subset"), col("layer"), col("cache_state"), col("phase"), col("rep")))
    end
    return keys
end

"""Append one result row. Keyword arguments override the defaults."""
function result_row!(ctx::Ctx; implementation, engine = "native", storage, workload, n_features = 0, gene_set = -1,
                     cell_subset = "", layer = "", n_cells_out = 0, phase = "fetch", rep = 0, wall_ms = NaN, cpu_ms = NaN,
                     open_ms = 0.0, alloc_mb = NaN, gc_ms = NaN, result_checksum = "", checksum_ok = "", status = "ok",
                     note = "", language = "julia", threads = Threads.nthreads())
    vals = Dict{String,Any}(
        "suite_version" => ctx.cfg["suite"]["version"], "run_id" => ctx.run_id,
        "timestamp" => Dates.format(now(UTC), "yyyy-mm-ddTHH:MM:SSZ"), "host" => ctx.host, "language" => language,
        "implementation" => implementation, "engine" => engine, "storage" => storage, "scale" => ctx.scale,
        "workload" => workload, "n_features" => n_features, "gene_set" => gene_set, "cell_subset" => cell_subset,
        "layer" => layer, "n_cells_out" => n_cells_out, "cache_state" => ctx.cache_state, "phase" => phase, "rep" => rep,
        "threads" => threads, "wall_ms" => _fmt(Float64(wall_ms)), "cpu_ms" => _fmt(Float64(cpu_ms)),
        "open_ms" => _fmt(Float64(open_ms)), "rss_mb" => _fmt(rss_mb()), "maxrss_mb" => _fmt(maxrss_mb()),
        "alloc_mb" => _fmt(Float64(alloc_mb)), "gc_ms" => _fmt(Float64(gc_ms)), "result_checksum" => result_checksum,
        "checksum_ok" => checksum_ok, "status" => status, "note" => note)
    println(ctx.io, join((_csv(vals[c]) for c in SCHEMA), ","))
    flush(ctx.io)
    push!(ctx.existing, row_key(implementation, workload, n_features, gene_set, cell_subset, layer, ctx.cache_state, phase, rep))
    return nothing
end

# ---------------------------------------------------------------------------
# time and memory

function cpu_time_ms()
    buf = zeros(UInt8, 256)
    @ccall getrusage(0::Cint, buf::Ptr{UInt8})::Cint
    if Sys.isapple()
        us = reinterpret(Int64, buf[1:8])[1] * 1_000_000 + reinterpret(Int32, buf[9:12])[1]
        ss = reinterpret(Int64, buf[17:24])[1] * 1_000_000 + reinterpret(Int32, buf[25:28])[1]
    else
        us = reinterpret(Int64, buf[1:8])[1] * 1_000_000 + reinterpret(Int64, buf[9:16])[1]
        ss = reinterpret(Int64, buf[17:24])[1] * 1_000_000 + reinterpret(Int64, buf[25:32])[1]
    end
    return (us + ss) / 1000
end

maxrss_mb() = Sys.maxrss() / 2^20

function rss_mb()
    try
        out = read(`ps -o rss= -p $(getpid())`, String)
        return parse(Float64, strip(out)) / 1024
    catch
        return NaN
    end
end

"""Run `f()` once; return (value, wall_ms, cpu_ms, alloc_mb, gc_ms)."""
function timed(f)
    c0 = cpu_time_ms()
    st = @timed f()
    c1 = cpu_time_ms()
    return st.value, st.time * 1000, c1 - c0, st.bytes / 2^20, st.gctime * 1000
end

# ---------------------------------------------------------------------------
# selections and checksums

function cells_for_subset(ctx::Ctx, all::Vector{String}, sub::AbstractString)
    sub == "all" && return nothing
    s = ctx.selection["cell_subsets"][sub]
    if startswith(sub, "contig")
        return all[(s[1] + 1):s[2]]
    end
    return all[Int.(s) .+ 1]
end

sha1_lines(v) = bytes2hex(sha1(join(string.(v), "\n")))

"""Float64 sum of the expression columns of a fetched table (all non-cell columns)."""
function feature_checksum(df::DataFrame)
    total = 0.0
    cc = SCUBA.cell_column(df)
    for n in names(df)
        n == cc && continue
        total += sum(Float64, df[!, n])
    end
    return total
end
feature_checksum(X::AbstractMatrix) = sum(Float64, X)
feature_checksum(X::SparseMatrixCSC) = sum(Float64, nonzeros(X))

metadata_checksum_full(df::DataFrame) = string(nrow(df), ":", repr(sum(Float64, df.nCount_RNA)))
function metadata_checksum(col::AbstractVector)
    if nonmissingtype(eltype(col)) <: Real && !(nonmissingtype(eltype(col)) <: Bool)
        return string(length(col), ":", repr(sum(Float64, skipmissing(col))))
    end
    return string(length(col), ":", sha1_lines(string.(col)))
end
reduction_checksum(df::DataFrame) = feature_checksum(df)
fetch_cells_checksum(ids) = string(length(ids), ":", sha1_lines(sort(String.(ids))))
unique_values_checksum(vals) = string(length(vals), ":", sha1_lines(sort(string.(vals))))

function checksum_ok(got, expected; rtol = 1e-5)
    expected === nothing && return ""
    if got isa Number && expected isa Number
        return string(isapprox(Float64(got), Float64(expected); rtol, atol = 1e-9))
    end
    if got isa AbstractString && expected isa AbstractString && occursin(':', got) && occursin(':', expected)
        a, b = split(got, ':', limit = 2), split(expected, ':', limit = 2)
        a[1] == b[1] || return "false"
        ga = tryparse(Float64, a[2]); gb = tryparse(Float64, b[2])
        (ga !== nothing && gb !== nothing) && return string(isapprox(ga, gb; rtol))
        return string(a[2] == b[2])
    end
    return string(got == expected)
end

# ---------------------------------------------------------------------------
# skip rules: "impl=<impl> scale=<scale> [n_features>=N] [workload=<w>]"

function skip_rule_matches(rules, impl, scale, n_features, workload)
    for r in rules
        parts = split(r)
        m = true
        for p in parts
            if startswith(p, "impl=")
                m &= (p[6:end] == impl)
            elseif startswith(p, "scale=")
                m &= (p[7:end] == scale)
            elseif startswith(p, "n_features>=")
                m &= (n_features >= parse(Int, p[13:end]))
            elseif startswith(p, "workload=")
                m &= (p[10:end] == workload)
            end
        end
        m && return r
    end
    return nothing
end

end # module
