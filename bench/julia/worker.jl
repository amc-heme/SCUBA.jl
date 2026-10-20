# Concurrency worker: runs a mixed query stream against one dataset for a fixed duration.
#
#   julia --project=bench -t 1 bench/julia/worker.jl --scale 1M --impl jl_h5ad_csc --worker 3 --duration 60
#         --out results/raw/<run_id>/conc/... --go <barrier file> [--seed 7] [--mode threads --sessions 8 --share shared|per_task]
#
# With --mode threads, this one process runs N tasks (Threads.@spawn) instead of being one of N processes.
# Each completed query appends a row (phase = concurrent) with wall_ms; a 5% sample is checksummed.
include(joinpath(@__DIR__, "Harness.jl"))
using .Harness
using SCUBA, DataFrames, OrderedCollections, JSON, Random, Dates
include(joinpath(@__DIR__, "common.jl"))

function wparse(args)
    o = Dict{String,Any}("config" => joinpath(Harness.bench_dir(), "config", "bench.toml"), "worker" => 1, "duration" => 60.0,
                         "seed" => 7, "mode" => "process", "n" => 1, "share" => "shared", "go" => "", "run_id" => "", "warmup" => 10.0)
    i = 1
    while i <= length(args)
        a = args[i]; v = args[i + 1]; i += 2
        if a == "--scale"; o["scale"] = v
        elseif a == "--impl"; o["impl"] = v
        elseif a == "--worker"; o["worker"] = parse(Int, v)
        elseif a == "--duration"; o["duration"] = parse(Float64, v)
        elseif a == "--warmup"; o["warmup"] = parse(Float64, v)
        elseif a == "--seed"; o["seed"] = parse(Int, v)
        elseif a == "--mode"; o["mode"] = v
        elseif a == "--sessions"; o["n"] = parse(Int, v)
        elseif a == "--share"; o["share"] = v
        elseif a == "--go"; o["go"] = v
        elseif a == "--run-id"; o["run_id"] = v
        elseif a == "--config"; o["config"] = v
        else error("unknown argument $a") end
    end
    return o
end

# Pick a workload from the configured mix
function pick(rng, ctx, mix_names, mix_w)
    r = rand(rng); acc = 0.0
    for (name, wgt) in zip(mix_names, mix_w)
        acc += wgt
        if r <= acc
            gs = ctx.selection["gene_sets"]
            if name == "feature1_all"; return Work("feature", "data", "1:$(rand(rng, 0:2))", 1, "all", "")
            elseif name == "feature10_rand10"; return Work("feature", "data", "10:$(rand(rng, 0:2))", 10, "rand10", "")
            elseif name == "meta_col"; return Work("metadata", "", "", 0, "all", rand(rng, ["cell_type", "nCount_RNA"]))
            elseif name == "reduction"; return Work("reduction", "", "", 0, "all", "X_umap|1-2")
            elseif name == "fetch_cells"; return Work("fetch_cells", "", "", 0, "", rand(rng, String.(ctx.cfg["workloads"]["fetch_cells"])))
            end
        end
    end
    return Work("feature", "data", "1:0", 1, "all", "")
end

function stream(ctx, src, worker_id, o, lock_io)
    rng = MersenneTwister(o["seed"] + worker_id)
    mix = ctx.cfg["concurrency"]["mix"]; names_ = collect(keys(mix)); ws = [mix[k] for k in names_]
    frac = ctx.cfg["concurrency"]["checksum_fraction"]
    all_ids = all_cells(src)
    # warm-up until the barrier file appears (or warmup seconds pass)
    t_warm = time()
    while isempty(o["go"]) ? (time() - t_warm < o["warmup"]) : !isfile(o["go"])
        run_work(ctx, src, pick(rng, ctx, names_, ws), all_ids)
        time() - t_warm > 300 && break
    end
    t_end = time() + o["duration"]
    n = 0
    while time() < t_end
        w = pick(rng, ctx, names_, ws)
        (chk, nout, exp), wall, cpu, alloc, gc = timed(() -> run_work(ctx, src, w, all_ids))
        ok = rand(rng) < frac ? checksum_ok(chk, exp) : ""
        lock(lock_io) do
            result_row!(ctx; implementation = o["impl"], storage = STORAGE[o["impl"]], workload = work_label(w),
                        n_features = w.n_features, gene_set = gene_set_index(w), cell_subset = work_subset_for_row(w),
                        layer = w.layer, n_cells_out = nout, phase = "concurrent", rep = worker_id, wall_ms = wall, cpu_ms = cpu,
                        alloc_mb = alloc, gc_ms = gc, result_checksum = chk isa Number ? repr(chk) : string(chk), checksum_ok = ok,
                        note = "mode=$(o["mode"]);share=$(o["share"]);n=$(o["n"])", threads = o["n"])
        end
        n += 1
    end
    return n
end

function main()
    o = wparse(ARGS)
    ctx = load_ctx(o["config"], o["scale"]; run_id = o["run_id"], cache_state = "warm")
    name = string("julia_conc_", o["impl"], "_", o["scale"], "_", o["mode"], "_n", o["n"], o["mode"] == "process" ? "_w$(o["worker"])" : "_$(o["share"])")
    Harness.open_output!(ctx, name)
    lk = ReentrantLock()
    if o["mode"] == "threads"
        src_shared = o["share"] == "shared" ? open_source(ctx, o["impl"]) : nothing
        tasks = [Threads.@spawn begin
                    src = src_shared === nothing ? open_source(ctx, o["impl"]) : src_shared
                    stream(ctx, src, i, o, lk)
                 end for i in 1:o["n"]]
        counts = fetch.(tasks)
        println("queries per worker: ", counts, " total ", sum(counts))
    else
        src = open_source(ctx, o["impl"])
        n = stream(ctx, src, o["worker"], o, lk)
        println("worker $(o["worker"]): $n queries")
    end
    close(ctx)
end

main()
