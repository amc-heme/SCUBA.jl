@testset "MetadataOverlay" begin
    src = synthetic_source()
    ov = MetadataOverlay(src)
    @test SCUBA.metadata_names(ov) == SCUBA.metadata_names(src)
    @test isequal(fetch_metadata(ov), fetch_metadata(src))
    @test isequal(fetch_data(ov, ["RNA_g1", "p1", "UMAP_1", "group", "score"]), fetch_data(src, ["RNA_g1", "p1", "UMAP_1", "group", "score"]))
    @test all_keys(ov) == all_keys(src)
    @test SCUBA.source_name(ov) == "MetadataOverlay(MemoryAnnData)"
    # replacing the table changes the overlay only
    SCUBA.set_metadata!(ov, DataFrame(new = collect(1:SCUBA.ncells(src)), grp = categorical(["x", "y", "x", "y", "x", "y"])))
    @test SCUBA.metadata_names(ov) == ["new", "grp"]
    @test fetch_data(ov, ["g1", "new"]).new == 1:6
    @test unique_values(ov, "grp") == ["x", "y"]
    @test "group" in SCUBA.metadata_names(src)
    @test_throws DimensionMismatch SCUBA.set_metadata!(ov, DataFrame(bad = 1:3))
    # an empty overlay hides every metadata column: only expression / reductions resolve
    ov2 = MetadataOverlay(src; obs = DataFrame())
    @test SCUBA.metadata_names(ov2) == String[]
    @test_throws VarsNotFound fetch_data(ov2, ["group"])
    @test names(fetch_data(ov2, ["g1", "UMAP_1"])) == ["cell", "g1", "UMAP_1"]
    # the subset of a conformance run over a file-backed inner source
    h5 = fixture_path("AML_h5ad.h5ad")
    if isfile(h5)
        inner = open_h5ad(h5)
        o = MetadataOverlay(inner)
        test_source_conformance(o; main_feature = "GAPDH", alt_modality = "protein", alt_feature = "CD117-AB",
                                reduction = "X_umap", categorical_var = "condensed_cell_type", numeric_var = "nCount_RNA")
        @test isopen(o)
        close(o)
        @test !isopen(o)
    end
end
