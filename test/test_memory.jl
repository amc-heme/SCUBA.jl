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

@testset "eltype option and dense-result guard" begin
    src = synthetic_source()
    df32 = fetch_data(src, ["g1_dummy_not_found", "RNA_g2", "AB_p1", "score"]; eltype = Float32, strict = false)
    @test eltype(df32.RNA_g2) == Float32 && eltype(df32.AB_p1) == Float32
    @test eltype(df32.score) == Float64                                   # metadata untouched
    @test df32.RNA_g2 == Float32.(fetch_data(src, ["RNA_g2"]).RNA_g2)
    @test eltype(fetch_feature(src, ["g2"]; eltype = Float32).g2) == Float32
    # the cell column shares the source's id vector when all cells are requested
    @test fetch_data(src, ["g2"]).cell === SCUBA.cell_ids(src)
    @test fetch_data(src, ["g2"]; cells = ["c2", "c1"]).cell == ["c2", "c1"]
    old = SCUBA.MAX_DENSE_BYTES[]
    SCUBA.MAX_DENSE_BYTES[] = 8 * 6 * 2 - 1       # two features × six cells of Float64 is just over
    try
        (_, w) = with_collected_warnings(() -> fetch_data(src, ["g2", "g3"]))
        @test any(x -> x.code == :large_result, w)
        (_, w2) = with_collected_warnings(() -> fetch_data(src, ["g2", "g3"]; eltype = Float32))
        @test !any(x -> x.code == :large_result, w2)
        SCUBA.DENSE_LIMIT_MODE[] = :error
        @test_throws ResultTooLarge fetch_data(src, ["g2", "g3"])
        @test fetch_data(src, ["g2"]).g2 == [i + 20 for i in 1:6]
    finally
        SCUBA.MAX_DENSE_BYTES[] = old
        SCUBA.DENSE_LIMIT_MODE[] = :warn
    end
    # column construction matches the matrix path on every synthetic layer
    for layer in ("data", "counts")
        cols = SCUBA.expression_columns(src, "RNA", layer, [2, 5, 1], [6, 1, 3], Float64)
        M = SCUBA.expression(src, "RNA", layer, [2, 5, 1], [6, 1, 3])
        @test reduce(hcat, cols) == M
    end
end
