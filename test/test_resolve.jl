@testset "resolve_vars rules" begin
    src = synthetic_source()

    @testset "keyed features keep their spelling; layer applies" begin
        df = fetch_data(src, ["RNA_g2", "AB_p1"])
        @test names(df) == ["cell", "RNA_g2", "AB_p1"]
        @test df.RNA_g2 == [i + 20 for i in 1:6]
        @test df.AB_p1 == [101 + i for i in 1:6]
        dfc = fetch_data(src, ["RNA_g2"]; layer = "counts")
        @test dfc.RNA_g2 == [2i for i in 1:6]
    end

    @testset "unkeyed main feature is bare; unkeyed alt feature is renamed with a warning" begin
        (df, w) = with_collected_warnings(() -> fetch_data(src, ["g2", "p1"]))
        @test names(df) == ["cell", "g2", "AB_p1"]
        @test length(w) == 1 && w[1].code == :alt_modality_renamed
    end

    @testset "metadata beats a same-named main feature, with a warning" begin
        (df, w) = with_collected_warnings(() -> fetch_data(src, ["g1"]))
        @test df.g1 == ["m$i" for i in 1:6]
        @test any(x -> x.code == :shadowed_feature, w)
        @test fetch_data(src, ["RNA_g1"]).RNA_g1 == [i + 10 for i in 1:6]
    end

    @testset "obs_ key addresses metadata and is stripped" begin
        df = fetch_data(src, ["obs_score", "obs_group"])
        @test names(df) == ["cell", "score", "group"]
        @test df.group isa CategoricalVector
        @test levels(df.group) == ["c", "b", "a"]
    end

    @testset "shared feature in several alternates: ambiguous only when not in main" begin
        # 'shared' is in RNA (main) -> main wins silently
        (df, w) = with_collected_warnings(() -> fetch_data(src, ["shared"]))
        @test names(df) == ["cell", "shared"]
        @test isempty(w)
    end

    @testset "key matching requires the underscore and is case-sensitive" begin
        @test_throws VarsNotFound fetch_data(src, ["RNAg1"])
        err = try fetch_data(src, ["rna_g1"]) catch e; e end
        @test err isa VarsNotFound
        @test err.suggestions["rna_g1"] == "RNA_g1"
    end

    @testset "reductions by dim; out of range is not found" begin
        df = fetch_data(src, ["UMAP_2", "PCA_3"])
        @test names(df) == ["cell", "UMAP_2", "PCA_3"]
        @test df.UMAP_2 == [i + 1.0 for i in 1:6]
        (df2, w) = with_collected_warnings(() -> fetch_data(src, ["UMAP_1", "UMAP_3"]))
        @test names(df2) == ["cell", "UMAP_1"]
        @test length(w) == 1 && w[1].code == :some_vars_not_found && w[1].vars == ["UMAP_3"]
        @test_throws VarsNotFound fetch_data(src, ["UMAP_1", "UMAP_3"]; strict = true)
    end

    @testset "duplicates collapse by target" begin
        (df, w) = with_collected_warnings(() -> fetch_data(src, ["g2", "g2", "RNA_g2", "p1", "AB_p1"]))
        @test names(df) == ["cell", "g2", "AB_p1"]
        @test count(x -> x.code == :duplicate_var, w) == 3
    end

    @testset "missing handling" begin
        @test_throws VarsNotFound fetch_data(src, ["nope"])
        (df, w) = with_collected_warnings(() -> fetch_data(src, ["g2", "nope", "nada"]))
        @test names(df) == ["cell", "g2"]
        @test w[end].code == :some_vars_not_found && w[end].vars == ["nope", "nada"]
        msg = sprint(showerror, VarsNotFound(["a"], 1, Dict{String,String}()))
        @test occursin("None of the requested variables were found", msg)
        msg2 = sprint(showerror, VarsNotFound(string.('a':'l'), 20, Dict{String,String}()))
        @test occursin("(10 out of 12 shown)", msg2)
    end

    @testset "column order follows input; cells subset and order" begin
        df = fetch_data(src, ["UMAP_1", "score", "g3", "AB_p2"]; cells = ["c4", "c1"])
        @test names(df) == ["cell", "UMAP_1", "score", "g3", "AB_p2"]
        @test df.cell == ["c4", "c1"]
        @test df.g3 == [34.0, 31.0]
        @test_throws CellsNotFound fetch_data(src, ["g3"]; cells = ["c1", "zz"])
    end

    @testset "layer substitution and missing layer" begin
        (df, w) = with_collected_warnings(() -> fetch_data(src, ["AB_p1"]; layer = "counts"))
        @test df.AB_p1 == [101 + i for i in 1:6]
        @test any(x -> x.code == :layer_substituted, w)
        @test_throws LayerNotFound fetch_data(src, ["g1", "RNA_g1"]; layer = "scaled")
    end

    @testset "nested keys: longest match with existence check" begin
        # modality named "X" plus a reduction named "X_umap", AnnData style
        cells = ["a", "b", "c"]
        X = sparse(Float64[1 2 3; 4 5 6; 7 8 9])
        s = MemoryAnnData(X; cells, features = ["umap_1", "g", "h"],
                          reductions = OrderedDict("X_umap" => Float64[10 20; 30 40; 50 60]))
        (df, w) = with_collected_warnings(() -> fetch_data(s, ["X_umap_1", "X_g"]))
        @test df.X_umap_1 == [10.0, 30.0, 50.0]      # reduction wins (longest key), but ambiguity is flagged
        @test any(x -> x.code == :ambiguous_var, w)
        @test df.X_g == [2.0, 5.0, 8.0]
        @test fetch_data(s, ["X_umap_2"]).X_umap_2 == [20.0, 40.0, 60.0]
    end

    @testset "metadata-only and reduction-only requests" begin
        df = fetch_data(src, ["group", "flag", "withmissing"])
        @test eltype(df.withmissing) == Union{Missing,Int}
        @test ismissing(df.withmissing[2])
        @test eltype(df.flag) == Bool
    end
end
