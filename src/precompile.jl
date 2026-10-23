# Precompile workload: exercise the public API on a tiny in-memory container (and on the bundled
# h5ad fixture when it is available locally) so first calls in a fresh session are fast. The work
# runs at precompile time only; nothing here touches the user's data.
using PrecompileTools

# Write a minimal anndata-layout HDF5 file so the backed reader's code paths precompile too.
function _write_tiny_h5ad(path, cells, genes, X::SparseMatrixCSC, obs::DataFrame, umap::Matrix{Float64})
    h5open(path, "w") do f
        attrs(f)["encoding-type"] = "anndata"; attrs(f)["encoding-version"] = "0.1.0"
        g = create_group(f, "X")
        attrs(g)["encoding-type"] = "csc_matrix"; attrs(g)["encoding-version"] = "0.1.0"; attrs(g)["shape"] = [size(X, 1), size(X, 2)]
        g["data"] = Float32.(nonzeros(X)); g["indices"] = Int32.(rowvals(X) .- 1); g["indptr"] = Int64.(X.colptr .- 1)
        for (name, index, cols) in (("obs", cells, obs), ("var", genes, DataFrame()))
            d = create_group(f, name)
            attrs(d)["encoding-type"] = "dataframe"; attrs(d)["encoding-version"] = "0.2.0"
            attrs(d)["_index"] = "_index"; attrs(d)["column-order"] = String[names(cols)...]
            d["_index"] = index
            for n in names(cols)
                v = cols[!, n]
                if v isa CategoricalVector
                    c = create_group(d, n); attrs(c)["encoding-type"] = "categorical"; attrs(c)["encoding-version"] = "0.2.0"; attrs(c)["ordered"] = false
                    c["codes"] = Int8[levelcode(x) - 1 for x in v]; c["categories"] = String.(levels(v))
                else
                    d[n] = v
                end
            end
        end
        m = create_group(f, "obsm"); attrs(m)["encoding-type"] = "dict"
        m["X_umap"] = permutedims(umap)
        attrs(m["X_umap"])["encoding-type"] = "array"
    end
end

"""
    _precompile_workload()

Run the API over tiny in-memory and on-disk containers. Used by PrecompileTools at precompile
time and by `tools/build_sysimage.jl` as the sysimage execution trace.
"""
function _precompile_workload()
    cells = ["c$i" for i in 1:8]
    rna = sparse(Float64[i % 3 == 0 ? 0 : i + j for i in 1:8, j in 1:5])
    ab = Float64[10 + i + j for i in 1:8, j in 1:3]
    obs = DataFrame(group = categorical(["a", "a", "b", "b", "c", "c", "a", "b"]), score = collect(1.0:8.0),
                    flag = [true, false, true, false, true, false, true, false], n = collect(1:8))
    reds = OrderedDict("UMAP" => Float64[i + 0.5j for i in 1:8, j in 1:2], "PC" => Float64[-i * j for i in 1:8, j in 1:3])
    begin
        src = MemoryAnnData(; cells, obs,
                            modalities = OrderedDict("RNA" => Modality(["g1", "g2", "g3", "g4", "g5"], OrderedDict("data" => rna, "counts" => rna)),
                                                     "AB" => Modality(["p1", "p2", "p3"], OrderedDict("data" => ab))),
                            reductions = reds, default_layer = "data")
        df = fetch_data(src, ["RNA_g1", "g2", "p1", "UMAP_1", "group", "score", "obs_flag", "nope"])
        fetch_data(src, ["g1"]; cells = ["c3", "c1"], layer = "counts")
        fetch_data(src, ["g1", "g2"]; eltype = Float32)
        fetch_feature(src, ["g1", "RNA_g3"])
        fetch_matrix(src, ["g1", "g2"]; sparse = true)
        fetch_matrix(src, ["g1", "g2"]; sparse = true, eltype = Float32)
        fetch_data(src, ["g1", "g2", "group", "UMAP_1"]; sparse = true)
        fetch_feature(src, ["g1"]; sparse = true, eltype = Float32)
        collect(eachchunk(src, ["g1", "UMAP_1", "group", "score"]; cells_per_chunk = 3))
        collect(eachchunk(src, ["g1"]; sparse = true, cells_per_chunk = 4))
        fetch_matrix(src, nothing)
        aggregate(src, Int32[1, 1, 2, 2, 3, 3, 1, 0], 3; layer = "counts")
        aggregate(src, [1, 2, 1], 2; cells = cells[1:3], features = ["g2", "g1"], fun = :mean, cells_per_chunk = 2)
        aggregate(src, Int32[1, 1, 2, 2, 1, 1, 2, 2], 2; modality = "AB")
        fetch_metadata(src); fetch_metadata(src, ["group", "score"]); fetch_metadata(src, "group")
        fetch_reduction(src); fetch_reduction(src, "PC"; dims = 1:3, cells = cells[1:3])
        fetch_cells(src, "group", ["a"]); unique_values(src, "group"); n_unique(src, "group")
        all_cells(src); meta_varnames(src); all_keys(src); features_in_modality(src); features_in_modality(src, "AB")
        reduction_dimnames(src, "UMAP", 1:2); default_layer(src); default_reduction(src)
        resolve_vars(src, ["g1", "UMAP_2"])
        with_collected_warnings(() -> fetch_data(src, ["g1", "missing_var"]))
        snapshot(src)
        SCUBA.cells(df); SCUBA.drop_cells(df); SCUBA.cell_column(df)
        tmp = joinpath(mktempdir(), "tiny.h5ad")
        try
            _write_tiny_h5ad(tmp, cells, ["g1", "g2", "g3", "g4", "g5"], rna, obs, reds["UMAP"])
            h5 = open_h5ad(tmp)
            fetch_data(h5, ["X_g1", "g2", "X_umap_1", "group", "score"])
            fetch_data(h5, ["g1"]; cells = cells[2:3], eltype = Float32)
            fetch_matrix(h5, ["g1", "g2"]; sparse = true); fetch_matrix(h5, ["g1"]; sparse = true, eltype = Float32)
            fetch_data(h5, ["X_g1", "group"]; sparse = true); fetch_metadata(h5); fetch_metadata(h5, "group")
            aggregate(h5, Int32[1, 2, 1, 2, 1, 2, 0, 0], 2); aggregate(h5, [1, 2], 2; cells = cells[1:2], features = ["g1"])
            fetch_reduction(h5); fetch_cells(h5, "group", ["a"]); unique_values(h5, "group")
            ov = MetadataOverlay(h5; obs = DataFrame())
            fetch_data(ov, ["X_g1", "g2"]; eltype = Float32); fetch_feature(ov, ["g1"]); fetch_reduction(ov)
            set_metadata!(ov, DataFrame(new = collect(1:8), grp = categorical(["a", "b", "a", "b", "a", "b", "a", "b"])))
            fetch_data(ov, ["g1", "new", "grp"]); fetch_metadata(ov); unique_values(ov, "grp")
            fetch_metadata(MetadataOverlay(h5), "group")
            snapshot(h5); close(h5)
        catch err
            @debug "precompile workload: h5ad part skipped" err
        end
        try
            fetch_data(src, ["definitely_missing"])
        catch e
            sprint(showerror, e)
        end
    end
end

@setup_workload begin
    @compile_workload begin
        _precompile_workload()
    end
end
