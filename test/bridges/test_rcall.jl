using RCall
@testset "RCall bridge (live Seurat / SingleCellExperiment)" begin
    have(pkg) = rcopy(Bool, rcall(:requireNamespace, pkg, quietly = true))
    h5 = fixture_path("AML_h5ad.h5ad")
    ref = open_h5ad(h5)
    vars_h5 = ["X_GAPDH", "protein_CD117-AB", "X_umap_1", "X_umap_2", "Batch", "nCount_RNA", "condensed_cell_type"]
    want = fetch_data(ref, vars_h5)
    # pseudobulk reference: counts summed per condensed_cell_type (sorted level order)
    ct = string.(fetch_data(ref, ["condensed_cell_type"]).condensed_cell_type)
    ct_levels = sort(unique(ct)); ng = length(ct_levels)
    ct_codes = Int32[findfirst(==(c), ct_levels) for c in ct]
    ind = zeros(Float64, length(ct), ng); for (i, g) in enumerate(ct_codes); ind[i, g] = 1.0; end
    fsh = features_in_modality(ref)
    want_counts = ind' * Matrix(fetch_matrix(ref, nothing; layer = "counts", sparse = true)[1])
    want_n = vec(sum(ind; dims = 1))
    align(V, fs) = V[:, [findfirst(==(f), fs) for f in fsh]]
    ScubaR = Base.get_extension(SCUBA, :SCUBARCallExt).ScubaR

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
            # aggregate: dgCMatrix layer scattered from R memory, in cell chunks bounded by CHUNK_BYTES
            V, n, fs, = aggregate(src, ct_codes, ng; layer = "counts")
            @test align(V, fs) ≈ want_counts && n == want_n
            @test SCUBA.aggregate_chunk_cells(src, "rna", "counts", length(fs)) isa Int
            old = SCUBA.CHUNK_BYTES[]; SCUBA.CHUNK_BYTES[] = 8 * length(fs) * 17     # 17 cells per chunk
            try
                @test align(aggregate(src, ct_codes, ng; layer = "counts")[1], fs) ≈ want_counts
            finally
                SCUBA.CHUNK_BYTES[] = old
            end
            # the R package's entry point: features × groups dense matrix
            pb = ScubaR.pseudobulk(reval("AML_Seurat"), ct_codes, ng, "RNA", "counts", nothing, nothing)
            @test align(permutedims(pb.values), pb.features) ≈ want_counts && pb.n == want_n
            pbm = ScubaR.pseudobulk(reval("AML_Seurat"), Int32[1, 2, 1], 2, nothing, "counts", all_cells(ref)[[3, 1, 2]], ["GAPDH"], "mean")
            g = fetch_data(ref, ["X_GAPDH"]; layer = "counts", cells = all_cells(ref)[[3, 1, 2]]).X_GAPDH
            @test vec(pbm.values) ≈ [(g[1] + g[3]) / 2, g[2]] && pbm.n == [2, 1] && pbm.features == ["GAPDH"]
            # The R helpers live in a private environment: clearing R's global environment (as
            # R CMD check does between examples, or `rm(list = ls())` in a session) must neither
            # break an existing source nor building a new one, and the helpers are not visible
            # in the user's workspace.
            helpers = ["scuba_slice", "scuba_seurat_layers", "scuba_seurat_layer_data", "scuba_bp_dir"]
            @test !any(h -> rcopy(Bool, rcall(:exists, h, envir = reval("globalenv()"), inherits = false)), helpers)
            reval("rm(list = ls(envir = globalenv(), all.names = TRUE), envir = globalenv())")
            @test fetch_data(src, ["rna_GAPDH"]; layer = "counts").rna_GAPDH == counts
            @test align(aggregate(src, ct_codes, ng; layer = "counts")[1], fs) ≈ want_counts
            fresh = rsource(reval("local({ load($(repr(rda))); AML_Seurat })"))
            @test fetch_data(fresh, ["rna_GAPDH"]; layer = "counts").rna_GAPDH == counts
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
            # aggregate streams the native BPCells layer (no chunking through R)
            @test SCUBA.aggregate_chunk_cells(src, "rna", "counts", 462) === nothing
            V, n, fs, = aggregate(src, ct_codes, ng; layer = "counts")
            @test align(V, fs) ≈ want_counts && n == want_n
            pb = ScubaR.pseudobulk(reval("AML_BPCells"), ct_codes, ng, "rna", "counts", nothing, nothing)
            @test align(permutedims(pb.values), pb.features) ≈ want_counts
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
            # aggregate over an HDF5Array-backed assay: slices are densified in R, so they are chunked
            @test SCUBA.aggregate_chunk_cells(src, "RNA", "counts", 462) isa Int
            V, n, fs, = aggregate(src, ct_codes, ng; layer = "counts")
            @test align(V, fs) ≈ want_counts && n == want_n
            pb = ScubaR.pseudobulk(reval("AML_SCE"), ct_codes, ng, "RNA", "counts", nothing, nothing)
            @test align(permutedims(pb.values), pb.features) ≈ want_counts
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
        # expression columns are R double vectors filled in place (zero-copy sink); `eltype` is accepted
        # but no longer changes the result
        rvec(c) = rcopy(Vector{Float64}, c)
        ref = open_h5ad(h5)
        res = ScubaR.fetch_data(h5, ["X_GAPDH", "protein_CD117-AB"], nothing, nothing, "Float32", false)
        @test res.names == ["X_GAPDH", "protein_CD117-AB"] && res.inputs == res.names
        @test res.columns[1] isa RObject && length(rvec(res.columns[1])) == 250
        GC.gc(true)                                   # the RObjects must survive a Julia GC before R reads them
        @test rvec(res.columns[1]) == fetch_data(ref, ["X_GAPDH"]).X_GAPDH
        @test rvec(res.columns[2]) == fetch_data(ref, ["protein_CD117-AB"])[!, "protein_CD117-AB"]
        res64 = ScubaR.fetch_data(h5, ["X_GAPDH"], nothing, ["487013_1", "39207_1"], "Float64", false)
        @test res64.cells == ["487013_1", "39207_1"] && res64.columns[1] isa RObject
        @test rvec(res64.columns[1]) == fetch_data(ref, ["X_GAPDH"]; cells = ["487013_1", "39207_1"]).X_GAPDH
        ff = ScubaR.fetch_feature(h5, ["GAPDH"], nothing, nothing, nothing, "Float64")
        @test ff.columns[1] isa RObject && rvec(ff.columns[1]) == fetch_data(ref, ["X_GAPDH"]).X_GAPDH
        # metadata / reduction columns still travel as Julia vectors
        red = ScubaR.fetch_reduction(h5, "X_umap", nothing, [1, 2])
        @test red.columns[1] isa Vector{Float64}
        close(ref)
        # R keeps the row names itself: `return_cells = false` omits the (potentially 1M-string) cell vector
        nocells = ScubaR.fetch_data(h5, ["X_GAPDH"], nothing, ["487013_1", "39207_1"], "Float64", false, false)
        @test nocells.cells === nothing && rvec(nocells.columns[1]) == rvec(res64.columns[1])
        @test ScubaR.fetch_feature(h5, ["GAPDH"], nothing, nothing, nothing, "Float64", false).cells === nothing
        @test ScubaR.fetch_reduction(h5, "X_umap", nothing, [1, 2], false).cells === nothing
        fm = ScubaR.fetch_matrix(h5, ["GAPDH"], nothing, nothing, nothing, true, "Float64", false)
        @test fm.cells === nothing && fm.dim == [1, 250]
        # pseudobulk over the file-backed path: R computes the codes, Julia returns features × groups
        ct = string.(fetch_data(open_h5ad(h5), ["condensed_cell_type"]).condensed_cell_type)
        lv = sort(unique(ct)); codes = Int32[findfirst(==(c), lv) for c in ct]
        ind = zeros(250, length(lv)); for (i, g) in enumerate(codes); ind[i, g] = 1.0; end
        r2 = open_h5ad(h5)
        pb = ScubaR.pseudobulk(h5, codes, length(lv), nothing, "counts", nothing, nothing)
        @test pb.values ≈ permutedims(ind' * Matrix(fetch_matrix(r2, nothing; layer = "counts", sparse = true)[1]))
        @test pb.n == vec(sum(ind; dims = 1)) && pb.features == features_in_modality(r2) && pb.warnings == String[]
        pbp = ScubaR.pseudobulk(h5, codes, length(lv), "protein", nothing, nothing, ["CD117-AB"], "mean")
        @test size(pbp.values) == (1, length(lv)) && pbp.features == ["CD117-AB"]
        close(r2)
        @test ScubaR.fetch_data(h5, ["Batch"], nothing, ["39207_1"], "Float64", true, false).cells === nothing
        # metadata is owned by R: a metadata name is not found in Julia; allow_empty yields an empty table
        empty = ScubaR.fetch_data(h5, ["Batch"], nothing, ["39207_1"], "Float64", true)
        @test empty.names == String[] && empty.cells == ["39207_1"] && occursin("Batch", only(empty.warnings))
        @test_throws SCUBA.VarsNotFound ScubaR.fetch_data(h5, ["Batch"], nothing, nothing, "Float64", false)
        # renamed alt-modality bare name keeps the requested spelling in `inputs`
        ren = ScubaR.fetch_data(h5, ["CD117-AB"], nothing, nothing, "Float64", false)
        @test ren.names == ["protein_CD117-AB"] && ren.inputs == ["CD117-AB"]
        ff = ScubaR.fetch_feature(h5, ["GAPDH"], nothing, nothing, nothing, "Float32")
        @test ff.names == ["GAPDH"] && ff.columns[1] isa RObject      # eltype accepted, result is an R double vector
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
        # expression! scatters straight from R memory (dgCMatrix slots / dense vector) into any sink
        n = length(SCUBA.all_cells(s)); gi = findfirst(==("GAPDH"), SCUBA.feature_names(s, "RNA"))
        cidx = collect(1:n); sub = [5, 1, 3]
        M = SCUBA.expression!(zeros(n, 1), s, "RNA", "data", [gi], cidx)
        @test vec(M) == d64[!, "rna_GAPDH"]
        cols = SCUBA.expression!(SCUBA._columns(Float32, 3, 1), s, "RNA", "data", [gi], sub)
        @test cols[1] == Float32.(d64[sub, "rna_GAPDH"])
        ai = findfirst(==("CD117-AB"), SCUBA.feature_names(s, "AB"))
        @test vec(SCUBA.expression!(zeros(n, 1), s, "AB", "data", [ai], cidx)) == d64[!, "ab_CD117-AB"]
        # the glue on a live object: R-vector columns, no cell vector
        g = ScubaR.fetch_data(reval("AML_Seurat"), ["rna_GAPDH"], nothing, nothing, "Float64", false, false)
        @test g.columns[1] isa RObject && g.cells === nothing
        GC.gc(true)
        @test rcopy(Vector{Float64}, g.columns[1]) == d64[!, "rna_GAPDH"]
    end
    if rcopy(Bool, rcall(:requireNamespace, "SingleCellExperiment", quietly = true)) && isfile(fixture_path("AML_SCE.rda"))
        reval("load($(repr(fixture_path("AML_SCE.rda"))))")
        sce = rsource(:AML_SCE)
        n = length(SCUBA.all_cells(sce)); gi = findfirst(==("GAPDH"), SCUBA.feature_names(sce, "RNA"))
        @test vec(SCUBA.expression!(zeros(n, 1), sce, "RNA", "logcounts", [gi], collect(1:n))) == fetch_data(sce, ["RNA_GAPDH"]).RNA_GAPDH
    end
end
