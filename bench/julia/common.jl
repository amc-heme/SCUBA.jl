# Source opening, workload enumeration and single-workload execution shared by run.jl and worker.jl.

const STORAGE = Dict("jl_h5ad_csc" => "h5ad_csc", "jl_h5ad_csr" => "h5ad_csr", "jl_bpcells" => "bpcells",
                     "jl_bpcells_rowmajor" => "bpcells_rowmajor", "jl_snapshot" => "memory")

# ---------------------------------------------------------------------------
# sources

function open_source(ctx, impl)
    d = ctx.data
    if impl == "jl_h5ad_csc"
        return open_h5ad(joinpath(d, "csc.h5ad"))
    elseif impl == "jl_h5ad_csr"
        return open_h5ad(joinpath(d, "csr.h5ad"))
    elseif impl == "jl_bpcells" || impl == "jl_bpcells_rowmajor"
        suffix = impl == "jl_bpcells" ? "" : "_rowmajor"
        h5 = open_h5ad(joinpath(d, "csc.h5ad"))
        obs = SCUBA.drop_cells(fetch_metadata(h5))
        reds = OrderedDict{String,Matrix{Float64}}()
        for r in SCUBA.reduction_keys(h5)
            reds[r] = Matrix(SCUBA.drop_cells(fetch_reduction(h5, r; dims = 1:SCUBA.reduction_ndims(h5, r))))
        end
        cells = all_cells(h5)
        close(h5)
        mod, bcells = SCUBA.bpcells_modality(Dict("counts" => joinpath(d, "bpcells", "RNA_counts" * suffix),
                                                  "data" => joinpath(d, "bpcells", "RNA_data" * suffix)))
        bcells == cells || error("BPCells cell names differ from the h5ad")
        return MemoryAnnData(; cells, obs, modalities = OrderedDict("RNA" => mod), reductions = reds, default_layer = "data")
    elseif impl == "jl_snapshot"
        h5 = open_h5ad(joinpath(d, "csc.h5ad"))
        s = snapshot(h5; layers = :all, chunk = 65536)
        close(h5)
        return s
    end
    error("unknown implementation $impl")
end

main_modality(src) = SCUBA.default_modality(src)
keyed(src, genes) = [string(main_modality(src), "_", g) for g in genes]

# ---------------------------------------------------------------------------
# workloads

struct Work
    workload::String
    layer::String
    gene_set::String        # "n:rep" or ""
    n_features::Int
    cell_subset::String
    extra::String           # reduction spec, metadata column, fetch_cells spec, ...
end

function workloads(ctx)
    wl = ctx.cfg["workloads"]
    W = Work[]
    for key in sort(collect(keys(ctx.selection["gene_sets"])); by = k -> (parse(Int, split(k, ":")[1]), parse(Int, split(k, ":")[2])))
        n, rep = parse.(Int, split(key, ":"))
        layers = copy(wl["layers_all_sets"])
        rep == 0 && append!(layers, wl["layers_first_set"])
        for layer in layers, sub in wl["cell_subsets"]
            push!(W, Work("feature", layer, key, n, sub, ""))
        end
        n == maximum(wl["n_features"]) && push!(W, Work("matrix_sparse", "data", key, n, "all", ""))
    end
    for m in wl["metadata"]
        push!(W, Work("metadata", "", "", 0, "all", m))
    end
    for r in wl["reductions"]
        push!(W, Work("reduction", "", "", 0, "all", string(r["name"], "|", r["dims"][1], "-", r["dims"][2])))
    end
    for f in wl["fetch_cells"]
        push!(W, Work("fetch_cells", "", "", 0, "", f))
    end
    for u in wl["unique_values"]
        push!(W, Work("unique_values", "", "", 0, "", u))
    end
    return W
end

function parse_single(spec)
    p = split(spec, "|")
    kind = p[1]
    if kind == "feature" || kind == "matrix_sparse"
        n = parse(Int, split(p[3], ":")[1])
        return Work(kind, p[2], p[3], n, p[4], "")
    elseif kind == "metadata"
        return Work("metadata", "", "", 0, "all", p[2])
    elseif kind == "reduction"
        return Work("reduction", "", "", 0, "all", p[2])
    else
        return Work(kind, "", "", 0, "", p[2])
    end
end

# run one workload once; returns (checksum, n_cells_out, expected)
function run_work(ctx, src, w::Work, all_ids)
    if w.workload == "feature"
        genes = String.(ctx.selection["gene_sets"][w.gene_set])
        cells = cells_for_subset(ctx, all_ids, w.cell_subset)
        layer = w.layer == "data" ? SCUBA.default_layer(src) : w.layer
        df = fetch_data(src, keyed(src, genes); layer, cells)
        e = get(ctx.expected["feature"], "$(w.layer)|$(w.gene_set)|$(w.cell_subset)", nothing)
        return feature_checksum(df), nrow(df), e === nothing ? nothing : e["sum"]
    elseif w.workload == "matrix_sparse"
        genes = String.(ctx.selection["gene_sets"][w.gene_set])
        X, c, g = fetch_matrix(src, genes; sparse = true)
        e = get(ctx.expected["feature"], "data|$(w.gene_set)|all", nothing)
        return feature_checksum(X), length(c), e === nothing ? nothing : e["sum"]
    elseif w.workload == "metadata"
        if w.extra == "full"
            df = fetch_metadata(src)
            return Harness.metadata_checksum_full(df), nrow(df), ctx.expected["metadata"]["full"]
        else
            col = fetch_metadata(src, w.extra)
            return metadata_checksum(col), length(col), ctx.expected["metadata"][w.extra]
        end
    elseif w.workload == "reduction"
        name, dims = split(w.extra, "|")
        lo, hi = parse.(Int, split(dims, "-"))
        df = fetch_reduction(src, name; dims = lo:hi)
        return reduction_checksum(df), nrow(df), ctx.expected["reduction"][w.extra]
    elseif w.workload == "fetch_cells"
        col, level = split(w.extra, ":")
        ids = fetch_cells(src, col, [level])
        return fetch_cells_checksum(ids), length(ids), ctx.expected["fetch_cells"][w.extra]
    elseif w.workload == "unique_values"
        v = unique_values(src, w.extra)
        return unique_values_checksum(v), length(v), ctx.expected["unique_values"][w.extra]
    end
    error("unknown workload $(w.workload)")
end

work_subset_for_row(w) = w.workload in ("feature", "matrix_sparse", "metadata", "reduction") ? w.cell_subset : ""

# Median warm wall time of the same cell at the next-smaller scale in this run, if recorded.
const SCALE_ORDER = ["10k", "100k", "1M", "3M"]
function previous_scale_medians(ctx, impl)
    i = findfirst(==(ctx.scale), SCALE_ORDER)
    (i === nothing || i == 1) && return Dict{String,Float64}(), 1.0
    prev = SCALE_ORDER[i - 1]
    path = joinpath(ctx.out_dir, string("julia_", impl, "_", prev, "_", ctx.cache_state, ".csv"))
    med = Dict{String,Vector{Float64}}()
    isfile(path) || return Dict{String,Float64}(), 1.0
    hdr = nothing
    for (k, line) in enumerate(eachline(path))
        if k == 1; hdr = split(line, ","); continue; end
        f = split(line, ","); length(f) >= length(hdr) || continue
        col(n) = f[findfirst(==(n), hdr)]
        (col("phase") == "fetch" && col("status") == "ok") || continue
        w = tryparse(Float64, col("wall_ms")); w === nothing && continue
        key = join((col("workload"), col("n_features"), col("cell_subset"), col("layer")), "|")
        push!(get!(med, key, Float64[]), w)
    end
    ratio = ctx.cfg["scales"][ctx.scale]["cells"] / ctx.cfg["scales"][prev]["cells"]
    return Dict(k => median(v) for (k, v) in med), ratio
end
gene_set_index(w) = isempty(w.gene_set) ? -1 : parse(Int, split(w.gene_set, ":")[2])
work_label(w) = w.workload in ("feature", "matrix_sparse") ? w.workload : string(w.workload, ":", w.extra)
note_for(w) = w.workload in ("feature", "matrix_sparse") ? "" : w.extra

