@testset "MuData (.h5mu) source" begin
    path = fixture_path("AML.h5mu")
    h5ad = fixture_path("AML_h5ad.h5ad")
    if isfile(path) && isfile(h5ad)
        mu = open_h5mu(path)
        @testset "schema" begin
            @test SCUBA.modality_keys(mu) == ["rna", "prot"]
            @test default_modality(mu) == "rna" && default_layer(mu) == "X"
            @test SCUBA.layer_names(mu, "rna") == ["X", "counts"] && SCUBA.layer_names(mu, "prot") == ["X"]
            @test length(all_cells(mu)) == 250
            @test Set(SCUBA.reduction_keys(mu)) == Set(["rna_X_pca", "rna_X_umap"])
            @test default_reduction(mu) == "rna_X_umap"
            @test SCUBA.reduction_dimnames(mu, "rna_X_umap", 1:2) == ["rna_X_umap_1", "rna_X_umap_2"]
            @test "condensed_cell_type" in meta_varnames(mu) && "rna_condensed_cell_type" in meta_varnames(mu)
            @test length(meta_varnames(mu)) == 56                       # 28 global + 28 rna-prefixed
            @test open_source(path) isa SCUBA.MuDataSource && open_source(h5ad) isa H5ADSource
            @test isopen(mu)
        end
        test_source_conformance(mu; main_feature = "GAPDH", alt_modality = "prot", alt_feature = "CD117-AB",
                                reduction = "rna_X_umap", categorical_var = "condensed_cell_type", numeric_var = "nCount_RNA")
        @testset "parity with the h5ad form" begin
            h = open_h5ad(h5ad)
            a = fetch_data(mu, ["rna_GAPDH", "prot_CD117-AB", "rna_X_umap_1", "Batch", "nCount_RNA"])
            b = fetch_data(h, ["X_GAPDH", "protein_CD117-AB", "X_umap_1", "Batch", "nCount_RNA"])
            @test a.cell == b.cell
            for (x, y) in zip(names(a)[2:end], names(b)[2:end])
                @test isequal(a[!, x], b[!, y])
            end
            @test fetch_data(mu, ["GAPDH"]; layer = "counts").GAPDH == fetch_data(h, ["GAPDH"]; layer = "counts").GAPDH
            @test isequal(fetch_data(mu, ["rna_Batch"]).rna_Batch, b.Batch)
            # bare reduction alias resolves to the only modality holding it; the column keeps the canonical name
            @test names(fetch_data(mu, ["X_umap_1"])) == ["cell", "rna_X_umap_1"]
            @test fetch_reduction(mu, "X_umap") == fetch_reduction(mu, "rna_X_umap")
            Xs, _, _ = fetch_matrix(mu, ["GAPDH"]; sparse = true)
            @test Vector(vec(Xs)) == b.X_GAPDH
            @test_throws SCUBA.MetadataNotFound fetch_metadata(mu, "prot_qc")    # prot obs is empty in this fixture
            close(h)
        end
        close(mu)
        @test !isopen(mu)
    end

    partial = fixture_path("AML_partial.h5mu")
    if isfile(partial) && isfile(h5ad)
        mu = open_h5mu(partial)
        full = open_h5ad(h5ad)
        cells = all_cells(mu)
        @testset "partial membership" begin
            @test length(cells) == 250
            absent_prot = findall(==(0), mu.modalities["prot"].obsmap)
            absent_rna = findall(==(0), mu.modalities["rna"].obsmap)
            @test length(absent_prot) == 50 && absent_rna == [248, 249, 250]
            d = fetch_data(mu, ["prot_CD117-AB", "rna_GAPDH"])
            @test findall(isnan, d[!, "prot_CD117-AB"]) == absent_prot
            @test findall(isnan, d[!, "rna_GAPDH"]) == absent_rna
            ref = fetch_data(full, ["protein_CD117-AB", "X_GAPDH"])
            present = setdiff(1:250, absent_prot)
            @test d[present, "prot_CD117-AB"] == ref[present, "protein_CD117-AB"]
            @test d[1:247, "rna_GAPDH"] == ref[1:247, "X_GAPDH"]
            # requested cell order mixing absent and present cells
            sub = cells[[absent_prot[1], 2, absent_prot[2], 7]]
            ds = fetch_data(mu, ["prot_CD117-AB"]; cells = sub)
            @test ds.cell == sub && isnan(ds[1, 2]) && isnan(ds[3, 2]) && !isnan(ds[2, 2])
            # sparse: absent cells have no structural entries
            Xs, _, _ = fetch_matrix(mu, ["CD117-AB"]; modality = "prot", sparse = true)
            @test all(iszero, Matrix(Xs)[absent_prot, 1])
            @test Matrix(Xs)[present, 1] == ref[present, "protein_CD117-AB"]
            # metadata: ambiguous column only keyed, modality-only column bare, absent cells missing
            @test "qc" ∉ meta_varnames(mu) && "rna_qc" in meta_varnames(mu) && "prot_qc" in meta_varnames(mu)
            @test "prot_only" in meta_varnames(mu)
            po = fetch_metadata(mu, "prot_only")
            @test findall(ismissing, po) == absent_prot
            @test SCUBA._is_categorical(po)
            @test findall(ismissing, fetch_metadata(mu, "rna_qc")) == absent_rna
            # reductions: modality-level with NaN rows, global one wins the bare alias
            r = fetch_reduction(mu, "prot_X_pca"; dims = 1:2)
            @test findall(isnan, r[!, 2]) == absent_prot
            @test "X_umap" in SCUBA.reduction_keys(mu)
            @test names(fetch_data(mu, ["X_umap_1"])) == ["cell", "X_umap_1"]
            @test fetch_data(mu, ["X_umap_1"]).X_umap_1 == fetch_data(full, ["X_umap_1"]).X_umap_1
            @test default_reduction(mu) == "X_umap"
            # chunks over a partial modality
            cat = reduce(vcat, collect(eachchunk(mu, ["prot_CD117-AB", "Batch"]; cells_per_chunk = 64)))
            @test isequal(cat[!, "prot_CD117-AB"], d[!, "prot_CD117-AB"])
            # overlay-backed metadata replacement
            o = MetadataOverlay(mu; obs = DataFrame(new = 1:250))
            @test fetch_data(o, ["rna_GAPDH", "new"]).new == 1:250
            @test_throws ArgumentError fetch_data(mu, ["rna_GAPDH"]; eltype = Int)
        end
        close(mu); close(full)
    end
end
