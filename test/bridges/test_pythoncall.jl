using PythonCall
@testset "PythonCall bridge (live anndata.AnnData)" begin
    ad = try
        pyimport("anndata")
    catch e
        @warn "python anndata not importable; skipping PythonCall bridge tests" exception = e
        nothing
    end
    h5 = fixture_path("AML_h5ad.h5ad")
    if ad !== nothing && isfile(h5)
        ref = open_h5ad(h5)
        adata = ad.read_h5ad(h5)
        src = pysource(adata)
        @test SCUBA.modality_keys(src) == ["X", "protein"]
        @test SCUBA.reduction_keys(src) == ["X_pca", "X_umap"]
        test_source_conformance(src; main_feature = "GAPDH", alt_modality = "protein", alt_feature = "CD117-AB",
                                reduction = "X_umap", categorical_var = "condensed_cell_type", numeric_var = "nCount_RNA")
        vars = ["X_GAPDH", "protein_CD117-AB", "X_umap_1", "X_umap_2", "Batch", "nCount_RNA", "ct"]
        a = fetch_data(src, vars); b = fetch_data(ref, vars)
        @test names(a) == names(b)
        for c in vars
            @test isequal(a[!, c], b[!, c]) || a[!, c] ≈ b[!, c]
        end
        @test levels(a.Batch) == levels(b.Batch)
        @test fetch_data(src, ["X_UNG"]; layer = "counts") == fetch_data(ref, ["X_UNG"]; layer = "counts")
        sub = all_cells(ref)[[9, 1, 250]]
        @test fetch_data(src, ["GAPDH", "X_umap_1"]; cells = sub) == fetch_data(ref, ["GAPDH", "X_umap_1"]; cells = sub)
        # aggregate: Python slices are dense, so the source is read in CHUNK_BYTES-bounded cell chunks
        ct = string.(fetch_data(ref, ["condensed_cell_type"]).condensed_cell_type)
        lv = sort(unique(ct)); codes = Int32[findfirst(==(c), lv) for c in ct]
        ind = zeros(250, length(lv)); for (i, g) in enumerate(codes); ind[i, g] = 1.0; end
        want_counts = ind' * Matrix(fetch_matrix(ref, nothing; layer = "counts", sparse = true)[1])
        @test SCUBA.aggregate_chunk_cells(src, "X", "counts", 462) isa Int
        V, n, = aggregate(src, codes, length(lv); layer = "counts")
        @test V ≈ want_counts && n == vec(sum(ind; dims = 1))
        @test aggregate(src, codes, length(lv); layer = "counts", cells_per_chunk = 9)[1] ≈ want_counts
        @test aggregate(src, codes, length(lv); modality = "protein")[1] ≈ ind' * Matrix(fetch_matrix(ref, nothing; modality = "protein")[1])

        @testset "backed mode" begin
            backed = ad.read_h5ad(h5; backed = "r")
            bs = pysource(backed)
            @test fetch_data(bs, ["X_GAPDH", "X_UNG"]) == fetch_data(ref, ["X_GAPDH", "X_UNG"])
            @test fetch_data(bs, ["X_GAPDH"]; cells = sub) == fetch_data(ref, ["X_GAPDH"]; cells = sub)
            @test aggregate(bs, codes, length(lv); layer = "counts")[1] ≈ want_counts
            backed.file.close()
        end

        @testset "set_metadata! writes through to Python" begin
            tbl = fetch_metadata(src)
            tbl.newcol = 1:250
            set_metadata!(src, tbl)
            @test pyconvert(Bool, "newcol" in adata.obs.columns)
            @test "newcol" in meta_varnames(src)
            @test fetch_data(src, ["newcol"]).newcol == 1:250
            @test string.(fetch_data(src, ["Batch"]).Batch) == string.(b.Batch)
        end
        close(ref)
    end
end

@testset "MuData parity with python mudata" begin
    mdmod = try
        pyimport("mudata")
    catch e
        @warn "python mudata not importable; skipping MuData parity tests" exception = e
        nothing
    end
    np = pyimport("numpy")
    partial = fixture_path("AML_partial.h5mu")
    if mdmod !== nothing && isfile(partial)
        m = mdmod.read_h5mu(partial)
        mu = open_h5mu(partial)
        obs_names = pyconvert(Vector{String}, m.obs_names.tolist())
        @test all_cells(mu) == obs_names
        # reindex a per-modality vector to the global cell order; absent cells -> NaN / missing
        function reindexed(mod, values, absent)
            names = pyconvert(Vector{String}, m.mod[mod].obs_names.tolist())
            pos = Dict(n => i for (i, n) in enumerate(names))
            [haskey(pos, c) ? values[pos[c]] : absent for c in obs_names]
        end
        function py_feature(mod, feat; layer = nothing)
            v = layer === nothing ? m.mod[mod].obs_vector(feat) : m.mod[mod].obs_vector(feat; layer)
            reindexed(mod, pyconvert(Vector{Float64}, np.asarray(v; dtype = "float64")), NaN)
        end
        d = fetch_data(mu, ["rna_GAPDH", "prot_CD117-AB", "rna_X_umap_1", "X_umap_2", "prot_only", "rna_qc"])
        @test isequal(d[!, "rna_GAPDH"], py_feature("rna", "GAPDH"))
        @test isequal(d[!, "prot_CD117-AB"], py_feature("prot", "CD117-AB"))
        @test isequal(fetch_data(mu, ["rna_GAPDH"]; layer = "counts")[!, "rna_GAPDH"], py_feature("rna", "GAPDH"; layer = "counts"))
        # modality-level reduction (NaN rows for absent cells) and the global obsm reduction
        umap = pyconvert(Matrix{Float64}, np.asarray(m.mod["rna"].obsm["X_umap"]; dtype = "float64"))
        @test isequal(d[!, "rna_X_umap_1"], reindexed("rna", umap[:, 1], NaN))
        gumap = pyconvert(Matrix{Float64}, np.asarray(m.obsm["X_umap"]; dtype = "float64"))
        @test d[!, "X_umap_2"] == gumap[:, 2]
        # per-modality obs columns: categorical and numeric, missing where the modality lacks the cell
        po = pyconvert(Vector{String}, m.mod["prot"].obs["prot_only"].astype("str").tolist())
        @test isequal(string.(coalesce.(d[!, "prot_only"], "missing")), reindexed("prot", po, "missing"))
        qc = pyconvert(Vector{Float64}, np.asarray(m.mod["rna"].obs["qc"]; dtype = "float64"))
        @test isequal(coalesce.(d[!, "rna_qc"], NaN), reindexed("rna", qc, NaN))
        # global obs columns come straight from the top-level obs
        @test string.(fetch_data(mu, ["Batch"]).Batch) == pyconvert(Vector{String}, m.obs["Batch"].astype("str").tolist())
        # obsmap agrees with mudata's
        for mod in ("rna", "prot")
            om = pyconvert(Vector{Int}, np.asarray(m.obsmap[mod]; dtype = "int64"))
            @test mu.modalities[mod].obsmap == om
        end
        close(mu)
    end
end
