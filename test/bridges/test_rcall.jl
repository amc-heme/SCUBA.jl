using RCall
@testset "RCall bridge (live Seurat / SingleCellExperiment)" begin
    have(pkg) = rcopy(Bool, rcall(:requireNamespace, pkg, quietly = true))
    h5 = fixture_path("AML_h5ad.h5ad")
    ref = open_h5ad(h5)
    vars_h5 = ["X_GAPDH", "protein_CD117-AB", "X_umap_1", "X_umap_2", "Batch", "nCount_RNA", "condensed_cell_type"]
    want = fetch_data(ref, vars_h5)

    function compare(df)
        @test df.cell == want.cell
        for j in 2:ncol(df)
            a = df[!, j]; b = want[!, j]
            if SCUBA._is_categorical(a)
                @test string.(a) == string.(b)
            else
                @test a ≈ b atol = 1e-4
            end
        end
    end

    if have("Seurat") && isfile(fixture_path("AML_Seurat.rda"))
        @testset "Seurat (in-memory v3 Assay)" begin
            rda = fixture_path("AML_Seurat.rda")
            reval("load($(repr(rda)))")
            src = rsource(:AML_Seurat)
            @test SCUBA.modality_keys(src) == ["rna", "ab"]
            @test SCUBA.reduction_keys(src) == ["PC", "UMAP"]
            @test default_reduction(src) == "UMAP"
            @test features_in_modality(src, "RNA") == features_in_modality(src, "rna")
            @test SCUBA.reduction_dimnames(src, "umap", [1, 2]) == ["UMAP_1", "UMAP_2"]
            test_source_conformance(src; main_feature = "GAPDH", alt_modality = "ab", alt_feature = "CD117-AB",
                                    reduction = "UMAP", categorical_var = "condensed_cell_type", numeric_var = "nCount_RNA")
            df = fetch_data(src, ["rna_GAPDH", "ab_CD117-AB", "UMAP_1", "UMAP_2", "Batch", "nCount_RNA", "condensed_cell_type"])
            compare(df)
            counts = fetch_data(src, ["rna_GAPDH"]; layer = "counts").rna_GAPDH
            @test counts == fetch_data(ref, ["X_GAPDH"]; layer = "counts").X_GAPDH
            @test_throws ReadOnlySource set_metadata!(src, fetch_metadata(src))
        end
    else
        @warn "SeuratObject not installed in R; skipping Seurat bridge tests"
    end

    if have("Seurat") && have("BPCells") && isfile(fixture_path("AML_BPCells.rda"))
        @testset "Seurat with BPCells layers (native fast path)" begin
            rda = fixture_path("AML_BPCells.rda")
            bpdir = fixture_path("AML_BPCells")
            reval("load($(repr(rda)))")
            # the stored matrices point at absolute paths from another machine: repoint them
            reval("""
            for (assay in c("RNA", "AB")) for (layer in c("counts", "data")) {
              d <- file.path($(repr(bpdir)), paste0(assay, "_", layer))
              AML_BPCells[[assay]]@layers[[layer]]@matrix@dir <- normalizePath(d)
            }
            """)
            src = rsource(:AML_BPCells)
            rna = src.modalities["rna"]
            @test Set(keys(rna.native)) == Set(["counts", "data"])
            df = fetch_data(src, ["rna_GAPDH", "ab_CD117-AB", "UMAP_1", "UMAP_2", "Batch", "nCount_RNA", "condensed_cell_type"])
            compare(df)
            @test fetch_data(src, ["rna_UNG"]; layer = "counts").rna_UNG == fetch_data(ref, ["X_UNG"]; layer = "counts").X_UNG
        end
    end

    if have("SingleCellExperiment") && have("HDF5Array") && isdir(fixture_path("AML_sce"))
        @testset "SingleCellExperiment (HDF5Array-backed)" begin
            d = fixture_path("AML_sce")
            reval("AML_SCE <- HDF5Array::loadHDF5SummarizedExperiment($(repr(d)))")
            src = rsource(:AML_SCE)
            @test SCUBA.modality_keys(src) == ["RNA", "AB"]
            @test SCUBA.reduction_keys(src) == ["PCA", "UMAP"]
            @test default_layer(src) == "logcounts"
            test_source_conformance(src; main_feature = "GAPDH", alt_modality = "AB", alt_feature = "CD117-AB",
                                    reduction = "UMAP", categorical_var = "condensed_cell_type", numeric_var = "nCount_RNA")
            df = fetch_data(src, ["RNA_GAPDH", "AB_CD117-AB", "UMAP_1", "UMAP_2", "Batch", "nCount_RNA", "condensed_cell_type"])
            compare(df)
            @test fetch_data(src, ["RNA_GAPDH"]; layer = "counts").RNA_GAPDH == fetch_data(ref, ["X_GAPDH"]; layer = "counts").X_GAPDH
            # single reduction dimension works (the R package mishandled this on SCE)
            @test names(fetch_data(src, ["UMAP_1"])) == ["cell", "UMAP_1"]
        end
    else
        @warn "SingleCellExperiment/HDF5Array not installed in R; skipping SCE bridge tests"
    end
    close(ref)
end

@testset "RCall bridge: default_reduction falls back to pca by name" begin
    if rcopy(Bool, rcall(:requireNamespace, "Seurat", quietly = true)) && isfile(fixture_path("AML_Seurat.rda"))
        reval("load($(repr(fixture_path("AML_Seurat.rda")))); sobj <- AML_Seurat; sobj[['umap']] <- NULL")
        s = rsource(:sobj)
        @test SCUBA.reduction_keys(s) == ["PC"]
        @test default_reduction(s) == "PC"
        @test SCUBA.reduction_search_names(s, "PC") == ["PC", "pca"]
    end
end

@testset "RCall bridge: file-backed glue entry points (ScubaR) and eltype" begin
    ScubaR = Base.get_extension(SCUBA, :SCUBARCallExt).ScubaR
    h5 = fixture_path("AML_h5ad.h5ad")
    if isfile(h5)
        ScubaR.clear_cache()
        info = ScubaR.source_open(h5)
        @test info.modalities[1] == "X" && "protein" in info.modalities
        @test length(info.features["X"]) == 462 && length(info.cells) == 250
        @test "X_umap" in info.reductions && info.default_reduction == "X_umap" && info.default_layer == "X"
        @test "condensed_cell_type" in info.obs.names && length(info.obs.cells) == 250
        @test info.uns_assays == String[]
        @test ScubaR.path_cache_size() == 1
        res = ScubaR.fetch_data(h5, ["X_GAPDH", "protein_CD117-AB"], nothing, nothing, "Float32", false)
        @test res.names == ["X_GAPDH", "protein_CD117-AB"] && res.inputs == res.names
        @test res.columns[1] isa Vector{Float32} && length(res.columns[1]) == 250
        res64 = ScubaR.fetch_data(h5, ["X_GAPDH"], nothing, ["487013_1", "39207_1"], "Float64", false)
        @test res64.cells == ["487013_1", "39207_1"] && res64.columns[1] isa Vector{Float64}
        # metadata is owned by R: a metadata name is not found in Julia; allow_empty yields an empty table
        empty = ScubaR.fetch_data(h5, ["Batch"], nothing, ["39207_1"], "Float64", true)
        @test empty.names == String[] && empty.cells == ["39207_1"] && occursin("Batch", only(empty.warnings))
        @test_throws SCUBA.VarsNotFound ScubaR.fetch_data(h5, ["Batch"], nothing, nothing, "Float64", false)
        # renamed alt-modality bare name keeps the requested spelling in `inputs`
        ren = ScubaR.fetch_data(h5, ["CD117-AB"], nothing, nothing, "Float64", false)
        @test ren.names == ["protein_CD117-AB"] && ren.inputs == ["CD117-AB"]
        ff = ScubaR.fetch_feature(h5, ["GAPDH"], nothing, nothing, nothing, "Float32")
        @test ff.names == ["GAPDH"] && ff.columns[1] isa Vector{Float32}
        @test_throws ArgumentError ScubaR.fetch_data(h5, ["X_GAPDH"], nothing, nothing, "Int64", false)
        @test ScubaR.path_cache_size() == 1
        ScubaR.source_close(h5)
        @test ScubaR.path_cache_size() == 0
    end
    if rcopy(Bool, rcall(:requireNamespace, "Seurat", quietly = true)) && isfile(fixture_path("AML_Seurat.rda"))
        reval("load($(repr(fixture_path("AML_Seurat.rda"))))")
        s = rsource(:AML_Seurat)
        d64 = fetch_data(s, ["rna_GAPDH", "ab_CD117-AB"])
        d32 = fetch_data(s, ["rna_GAPDH", "ab_CD117-AB"]; eltype = Float32)
        @test eltype(d32[!, "rna_GAPDH"]) == Float32 && eltype(d32[!, "ab_CD117-AB"]) == Float32
        @test d32[!, "rna_GAPDH"] ≈ Float32.(d64[!, "rna_GAPDH"])
        @test d32[!, "ab_CD117-AB"] ≈ Float32.(d64[!, "ab_CD117-AB"])
    end
end
