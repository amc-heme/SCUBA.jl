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

    if have("SeuratObject") && isfile(fixture_path("AML_Seurat.rda"))
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

    if have("SeuratObject") && have("BPCells") && isfile(fixture_path("AML_BPCells.rda"))
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
    if rcopy(Bool, rcall(:requireNamespace, "SeuratObject", quietly = true)) && isfile(fixture_path("AML_Seurat.rda"))
        reval("load($(repr(fixture_path("AML_Seurat.rda")))); sobj <- AML_Seurat; sobj[['umap']] <- NULL")
        s = rsource(:sobj)
        @test SCUBA.reduction_keys(s) == ["PC"]
        @test default_reduction(s) == "PC"
        @test SCUBA.reduction_search_names(s, "PC") == ["PC", "pca"]
    end
end
