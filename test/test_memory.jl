@testset "MemoryAnnData" begin
    src = synthetic_source()
    test_source_conformance(src; main_feature = "g2", alt_modality = "AB", alt_feature = "p1", reduction = "UMAP",
                            categorical_var = "group", numeric_var = "score")
    @testset "construction checks" begin
        @test_throws DimensionMismatch MemoryAnnData(zeros(3, 2); cells = ["a", "b"], features = ["x", "y"])
        @test_throws DimensionMismatch MemoryAnnData(zeros(2, 2); cells = ["a", "b"], features = ["x", "y"],
                                                     reductions = OrderedDict("U" => zeros(3, 2)))
    end
    @testset "set_metadata!" begin
        tbl = fetch_metadata(src)
        tbl.newcol = 1:6
        set_metadata!(src, tbl)
        @test "newcol" in meta_varnames(src)
        @test !("cell" in meta_varnames(src))
        @test fetch_data(src, ["newcol"]).newcol == 1:6
        @test_throws DimensionMismatch set_metadata!(src, DataFrame(a = 1:2))
    end
    @testset "dense layers and LazyMatrix" begin
        struct DoublingLazy <: SCUBA.LazyMatrix
            A::Matrix{Float64}
        end
        Base.size(m::DoublingLazy) = size(m.A)
        SCUBA.pull(m::DoublingLazy, fidx, cidx) = 2 .* m.A[cidx, fidx]
        A = Float64[i * 10 + j for i in 1:4, j in 1:3]
        s = MemoryAnnData(; cells = ["a", "b", "c", "d"],
                          modalities = OrderedDict("M" => Modality(["f1", "f2", "f3"], OrderedDict("X" => DoublingLazy(A), "dense" => A))))
        @test fetch_data(s, ["f2"]).f2 == 2 .* A[:, 2]
        @test fetch_data(s, ["f2"]; layer = "dense").f2 == A[:, 2]
        @test Matrix(fetch_matrix(s, ["f3", "f1"]; sparse = true)[1]) == 2 .* A[:, [3, 1]]
    end
    @testset "default_reduction preference" begin
        mk(reds) = MemoryAnnData(zeros(2, 1); cells = ["a", "b"], features = ["g"],
                                 reductions = OrderedDict(r => zeros(2, 2) for r in reds))
        @test default_reduction(mk(["pca", "tsne", "umap"])) == "umap"
        @test default_reduction(mk(["PCA", "TSNE"])) == "TSNE"
        @test default_reduction(mk(["X_pca"])) == "X_pca"
        @test_throws NoDefaultReduction default_reduction(mk(String[]))
        @test_throws NoDefaultReduction default_reduction(mk(["harmony"]))
    end
    @testset "cell column collision" begin
        s = MemoryAnnData(zeros(2, 1); cells = ["a", "b"], features = ["cell"], obs = DataFrame(cell_id = [1, 2]))
        df = fetch_data(s, ["cell", "cell_id"])
        @test cell_column(df) == "scuba_cell"
        @test cells(df) == ["a", "b"]
        @test names(drop_cells(df)) == ["cell", "cell_id"]
    end
end
