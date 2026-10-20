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

        @testset "backed mode" begin
            backed = ad.read_h5ad(h5; backed = "r")
            bs = pysource(backed)
            @test fetch_data(bs, ["X_GAPDH", "X_UNG"]) == fetch_data(ref, ["X_GAPDH", "X_UNG"])
            @test fetch_data(bs, ["X_GAPDH"]; cells = sub) == fetch_data(ref, ["X_GAPDH"]; cells = sub)
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
