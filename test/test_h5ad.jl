@testset "H5ADSource (AML fixture)" begin
    path = fixture_path("AML_h5ad.h5ad")
    if !isfile(path)
        @warn "AML_h5ad.h5ad fixture not found; skipping" path
    else
        src = open_h5ad(path)
        @test SCUBA.modality_keys(src) == ["X", "protein"]
        @test SCUBA.layer_names(src, "X") == ["X", "counts"]
        @test SCUBA.reduction_keys(src) == ["X_pca", "X_umap"]
        @test length(all_cells(src)) == 250
        @test length(features_in_modality(src)) == 462
        @test length(features_in_modality(src, "protein")) == 197
        test_source_conformance(src; main_feature = "GAPDH", alt_modality = "protein", alt_feature = "CD117-AB",
                                reduction = "X_umap", categorical_var = "condensed_cell_type", numeric_var = "nCount_RNA")
        # AnnData-specific: X_ key kept, obs_ stripped, counts layer differs from X
        df = fetch_data(src, ["X_GAPDH", "obs_Batch", "X_umap_2"])
        @test names(df) == ["cell", "X_GAPDH", "Batch", "X_umap_2"]
        @test levels(df.Batch) == ["BM_200AB", "PBMC_200AB"]
        counts = fetch_data(src, ["X_GAPDH"]; layer = "counts").X_GAPDH
        @test all(isinteger, counts)
        @test counts != df.X_GAPDH
        @test expm1.(df.X_GAPDH) ≈ counts ./ fetch_data(src, ["nCount_RNA"]).nCount_RNA .* 1e4 rtol = 1e-6
        # CSC vs sparse extraction agree
        X, _, _ = fetch_matrix(src, ["GAPDH", "UNG", "MEIS1"]; sparse = true)
        @test Matrix(X) == Matrix(fetch_matrix(src, ["GAPDH", "UNG", "MEIS1"])[1])
        # full-cells and partial-cells reads agree
        sub = all_cells(src)[200:-7:5]
        @test fetch_data(src, ["X_UNG"]; cells = sub).X_UNG == fetch_data(src, ["X_UNG"]).X_UNG[200:-7:5]
        @test fetch_cells(src, "Batch", "PBMC_200AB") == all_cells(src)[fetch_data(src, ["Batch"]).Batch .== "PBMC_200AB"]
        close(src)
        @test !isopen(src)
    end
end

@testset "H5AD variants (csr, dense, edge cases)" begin
    base = fixture_path("AML_h5ad.h5ad")
    csr_p = fixture_path("AML_h5ad_csr.h5ad")
    dense_p = fixture_path("AML_h5ad_dense.h5ad")
    edge_p = fixture_path("AML_h5ad_edge.h5ad")
    if !(isfile(base) && isfile(csr_p) && isfile(dense_p) && isfile(edge_p))
        @warn "derived h5ad fixtures not found (run fixtures/make_fixtures.py); skipping"
    else
        ref = open_h5ad(base)
        vars = ["X_GAPDH", "X_UNG", "protein_CD117-AB", "X_umap_1", "Batch"]
        want = fetch_data(ref, vars)
        want_counts = fetch_data(ref, ["X_GAPDH", "X_UNG"]; layer = "counts")
        sub = all_cells(ref)[[7, 250, 1, 128]]
        want_sub = fetch_data(ref, vars; cells = sub)
        for (label, p) in (("csr", csr_p), ("dense", dense_p))
            @testset "$label" begin
                s = open_h5ad(p)
                kind = SCUBA._h5ref(s, "X", "X").kind
                @test kind == (label == "csr" ? :csr : :dense)
                got = fetch_data(s, vars)
                @test names(got) == names(want)
                for c in vars
                    @test isequal(got[!, c], want[!, c]) || got[!, c] ≈ want[!, c]
                end
                @test fetch_data(s, ["X_GAPDH", "X_UNG"]; layer = "counts") == want_counts
                @test fetch_data(s, vars; cells = sub)[!, "X_GAPDH"] == want_sub[!, "X_GAPDH"]
                # sparse extraction path
                Xs, _, _ = fetch_matrix(s, ["GAPDH", "UNG"]; sparse = true)
                @test Matrix(Xs) == Matrix(want[!, ["X_GAPDH", "X_UNG"]])
                # a very small read limit forces multi-block reads
                old = SCUBA.MAX_READ_ELEMENTS[]
                SCUBA.MAX_READ_ELEMENTS[] = 100
                try
                    @test fetch_data(s, ["X_GAPDH", "X_UNG", "X_MEIS1"]).X_UNG == want.X_UNG
                    @test fetch_data(s, ["X_GAPDH"]; cells = sub).X_GAPDH == want_sub.X_GAPDH
                finally
                    SCUBA.MAX_READ_ELEMENTS[] = old
                end
                close(s)
            end
        end
        @testset "edge-case metadata encodings" begin
            s = open_h5ad(edge_p)
            df = fetch_data(s, ["unused_level", "with_na", "nullable_int", "flag"])
            @test df.unused_level isa CategoricalVector
            @test levels(df.unused_level) == ["BM_200AB", "PBMC_200AB", "UNUSED"]
            @test string.(df.unused_level) == string.(fetch_data(ref, ["Batch"]).Batch)
            @test eltype(df.with_na) == Float64 && isnan(df.with_na[1]) && !isnan(df.with_na[2])
            @test eltype(df.nullable_int) == Union{Missing,Int64}
            @test ismissing(df.nullable_int[1]) && df.nullable_int[2] == 1
            @test eltype(df.flag) == Bool && df.flag[1] == true && df.flag[2] == false
            @test fetch_cells(s, "flag", [true]) == all_cells(s)[1:2:end]
            @test n_unique(s, "unused_level") == 2
            close(s)
        end
        close(ref)
    end
end

@testset "H5AD column construction matches the matrix path" begin
    path = fixture_path("AML_h5ad.h5ad")
    if isfile(path)
        src = open_h5ad(path)
        fidx = [10, 3, 400]; cidx = [250, 1, 77, 78]
        for layer in ("X", "counts")
            cols = SCUBA.expression_columns(src, "X", layer, fidx, cidx, Float64)
            @test reduce(hcat, cols) == SCUBA.expression(src, "X", layer, fidx, cidx)
            @test eltype(SCUBA.expression_columns(src, "X", layer, fidx, cidx, Float32)[1]) == Float32
        end
        colsp = SCUBA.expression_columns(src, "protein", "X", [1, 5], cidx, Float64)
        @test reduce(hcat, colsp) == SCUBA.expression(src, "protein", "X", [1, 5], cidx)
        df = fetch_data(src, ["X_GAPDH", "protein_CD117-AB"]; eltype = Float32)
        @test eltype(df.X_GAPDH) == Float32 && eltype(df[!, "protein_CD117-AB"]) == Float32
        close(src)
    end
end
