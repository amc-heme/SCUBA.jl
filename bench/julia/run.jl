# Julia benchmark harness for SCUBA.jl native readers.
#
#   julia --project=bench -t 1 bench/julia/run.jl --scale 10k [--config bench/config/bench.toml]
#        [--cache warm|cold] [--only jl_h5ad_csc,jl_bpcells] [--include-optional] [--resume]
#        [--run-id ID] [--single "feature|data|10:0|all" --rep N] [--sweep max_read_elements=1048576,8388608]
#        [--start-time <unix seconds>]   (passed by the Makefile to record process startup)
#
# One row per repetition is appended to results/raw/<run_id>/julia_<impl>_<scale>_<cache>.csv.
const T_PROCESS_START = time()
using Pkg
include(joinpath(@__DIR__, "Harness.jl"))
using .Harness
using SCUBA, DataFrames, CategoricalArrays, SparseArrays, OrderedCollections, JSON, Statistics
const T_LOADED = time()

# ---------------------------------------------------------------------------
# arguments

function parse_args(args)
    opts = Dict{String,Any}("config" => joinpath(Harness.bench_dir(), "config", "bench.toml"), "cache" => "warm",
                            "only" => String[], "include_optional" => false, "resume" => false, "run_id" => "",
                            "single" => "", "rep" => 1, "sweep" => "", "start_time" => NaN)
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--include-optional"; opts["include_optional"] = true
        elseif a == "--resume"; opts["resume"] = true
        else
            v = args[i + 1]; i += 1
            if a == "--scale"; opts["scale"] = v
            elseif a == "--config"; opts["config"] = v
            elseif a == "--cache"; opts["cache"] = v
            elseif a == "--only"; opts["only"] = split(v, ",")
            elseif a == "--run-id"; opts["run_id"] = v
            elseif a == "--single"; opts["single"] = v
            elseif a == "--rep"; opts["rep"] = parse(Int, v)
            elseif a == "--sweep"; opts["sweep"] = v
            elseif a == "--start-time"; opts["start_time"] = parse(Float64, v)
            else error("unknown argument $a") end
        end
        i += 1
    end
    haskey(opts, "scale") || error("--scale is required")
    return opts
end

include(joinpath(@__DIR__, "common.jl"))

function run_impl(ctx, impl, opts)
    storage = STORAGE[impl]
    tun = "max_read_elements=$(SCUBA.MAX_READ_ELEMENTS[]);bp_chunk_cache=$(SCUBA.BP_CHUNK_CACHE_MAX[])"
    reps = ctx.cfg["reps"]
    warm = reps["warm"]
    timeout_ms = reps["timeout_s"] * 1000
    rules = ctx.cfg["skip"]["rules"]

    # open
    src = nothing
    if !(opts["resume"] && row_key(impl, "open", 0, -1, "", "", ctx.cache_state, "open", 0) in ctx.existing && isempty(opts["single"]))
        src, wall, cpu, alloc, gc = timed(() -> open_source(ctx, impl))
        result_row!(ctx; implementation = impl, storage, workload = "open", phase = "open", rep = 0, wall_ms = wall,
                    cpu_ms = cpu, open_ms = wall, alloc_mb = alloc, gc_ms = gc, note = tun)
    else
        src = open_source(ctx, impl)
    end
    all_ids = all_cells(src)

    works = isempty(opts["single"]) ? workloads(ctx) : [parse_single(opts["single"])]
    prev_med, ratio = previous_scale_medians(ctx, impl)
    factor = get(ctx.cfg["skip"], "extrapolation_factor", 2.0)
    slow_ms = get(reps, "reduce_reps_above_s", 30) * 1000
    seen_family = Set{String}()
    for w in works
        rule = skip_rule_matches(rules, impl, ctx.scale, w.n_features, w.workload)
        if rule !== nothing
            result_row!(ctx; implementation = impl, storage, workload = work_label(w), n_features = w.n_features,
                        gene_set = gene_set_index(w), cell_subset = work_subset_for_row(w), layer = w.layer,
                        status = "skipped", note = "rule: $rule")
            continue
        end
        pkey = join((work_label(w), w.n_features, work_subset_for_row(w), w.layer), "|")
        if haskey(prev_med, pkey) && prev_med[pkey] * ratio * factor > timeout_ms
            result_row!(ctx; implementation = impl, storage, workload = work_label(w), n_features = w.n_features,
                        gene_set = gene_set_index(w), cell_subset = work_subset_for_row(w), layer = w.layer,
                        status = "skipped_extrapolated",
                        note = "predicted $(round(prev_med[pkey] * ratio / 1000; digits = 1)) s from previous scale")
            continue
        end
        family = string(w.workload, ":", w.layer)
        nrep = isempty(opts["single"]) ? warm : 1
        rep_start = isempty(opts["single"]) ? 0 : opts["rep"]
        budget_exceeded = false
        rep = rep_start - 1
        while (rep += 1) <= rep_start + nrep - (isempty(opts["single"]) ? 0 : 1)
            phase = rep == 0 ? (family in seen_family ? "warmup" : "first_call") : "fetch"
            isempty(opts["single"]) && rep == 0 && push!(seen_family, family)
            key = row_key(impl, work_label(w), w.n_features, gene_set_index(w), work_subset_for_row(w), w.layer, ctx.cache_state, phase, rep)
            opts["resume"] && key in ctx.existing && continue
            if budget_exceeded
                result_row!(ctx; implementation = impl, storage, workload = work_label(w), n_features = w.n_features,
                            gene_set = gene_set_index(w), cell_subset = work_subset_for_row(w), layer = w.layer, phase, rep,
                            status = "skipped_extrapolated", note = "first call exceeded timeout budget")
                continue
            end
            GC.gc(false)
            local chk, nout, exp, wall, cpu, alloc, gc, status, note
            try
                (chk, nout, exp), wall, cpu, alloc, gc = timed(() -> run_work(ctx, src, w, all_ids))
                status = "ok"; note = note_for(w)
            catch err
                chk = ""; nout = 0; exp = nothing; wall = NaN; cpu = NaN; alloc = NaN; gc = NaN
                status = "error"; note = first(sprint(showerror, err), 300)
            end
            ok = status == "ok" ? checksum_ok(chk, exp) : ""
            result_row!(ctx; implementation = impl, storage, workload = work_label(w), n_features = w.n_features,
                        gene_set = gene_set_index(w), cell_subset = work_subset_for_row(w), layer = w.layer,
                        n_cells_out = nout, phase, rep, wall_ms = wall, cpu_ms = cpu, alloc_mb = alloc, gc_ms = gc,
                        result_checksum = chk isa Number ? repr(chk) : string(chk), checksum_ok = ok, status, note)
            if status == "ok" && ok == "false"
                @warn "checksum mismatch" impl workload = work_label(w) gene_set = w.gene_set subset = w.cell_subset got = chk expected = exp
            end
            if status == "ok" && wall > timeout_ms
                budget_exceeded = true
            elseif status == "ok" && rep == 0 && wall > slow_ms
                nrep = 1                         # slow cell: one recorded repetition is enough
            end
        end
    end
    src isa SCUBA.H5ADSource && close(src)
    return nothing
end

function main()
    opts = parse_args(ARGS)
    ctx = load_ctx(opts["config"], opts["scale"]; run_id = opts["run_id"], cache_state = opts["cache"])
    impls = isempty(opts["only"]) ? copy(ctx.cfg["implementations"]["julia"]) : String.(opts["only"])
    if isempty(opts["only"]) && !opts["include_optional"]
        filter!(i -> !(i in ctx.cfg["implementations"]["optional"]), impls)
    end
    sweeps = isempty(opts["sweep"]) ? [nothing] : parse.(Int, split(split(opts["sweep"], "=")[2], ","))
    for sweep in sweeps
        sweep === nothing || (SCUBA.MAX_READ_ELEMENTS[] = sweep)
        for impl in impls
            name = string("julia_", impl, "_", ctx.scale, "_", ctx.cache_state, sweep === nothing ? "" : "_mre$sweep")
            Harness.open_output!(ctx, name)
            # startup rows (once per process, attached to the first implementation)
            if impl == impls[1] && sweep === nothing && isempty(opts["single"])
                t0 = isnan(opts["start_time"]) ? T_PROCESS_START : opts["start_time"]
                result_row!(ctx; implementation = impl, storage = STORAGE[impl], workload = "startup", phase = "startup",
                            rep = 0, wall_ms = (T_LOADED - t0) * 1000, note = "process start -> using SCUBA (incl. Julia startup when --start-time given)")
            end
            @info "running" impl scale = ctx.scale cache = ctx.cache_state threads = Threads.nthreads() sweep
            run_impl(ctx, impl, opts)
            close(ctx)
        end
    end
    println("results in ", ctx.out_dir)
end

main()
