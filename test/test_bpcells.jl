@testset "BPCells directory reader" begin
    dir = fixture_path("AML_BPCells")
    h5 = fixture_path("AML_h5ad.h5ad")
    if !(isdir(dir) && isfile(h5))
        @warn "BPCells fixtures not found; skipping"
    else
        @testset "bp128 decoder round trips" begin
            # pack with a simple reference packer and check unpack
            function refpack(vals::Vector{UInt32}, bits)
                words = zeros(UInt32, 4 * bits)
                for lane in 0:3, i in 0:31
                    v = vals[4i + lane + 1]
                    bitpos = i * bits
                    for b in 0:(bits - 1)
                        bit = (v >> b) & 1
                        pos = bitpos + b
                        words[4 * (pos ÷ 32) + lane + 1] |= bit << (pos % 32)
                    end
                end
                words
            end
            for bits in (1, 3, 7, 8, 13, 17, 31)
                vals = UInt32.(rand(0:(2^bits - 1), 128))
                out = SCUBA.bp128_unpack!(zeros(UInt32, 128), refpack(vals, bits), 0, bits)
                @test out == vals
            end
            @test SCUBA.bp128_unpack!(ones(UInt32, 128), UInt32[], 0, 0) == zeros(UInt32, 128)
            w = UInt32.(1:128)
            @test SCUBA.bp128_unpack!(zeros(UInt32, 128), w, 0, 32) == w
        end

        rna_data = open_bpcells(joinpath(dir, "RNA_data"))
        rna_counts = open_bpcells(joinpath(dir, "RNA_counts"))
        ab_data = open_bpcells(joinpath(dir, "AB_data"))
        @test size(rna_data) == (250, 462)
        @test rna_data.version == "packed-double-matrix-v2"
        @test SCUBA.nnz_entries(rna_data) == 33374
        @test length(SCUBA.row_names(rna_data)) == 462 && length(SCUBA.col_names(rna_data)) == 250

        src = open_h5ad(h5)
        @test SCUBA.col_names(rna_data) == all_cells(src)
        @test SCUBA.row_names(rna_data) == features_in_modality(src)
        # the BPCells RNA 'data' layer is the same log-normalised matrix as the h5ad X
        X_bp = Matrix(SparseArrays.sparse(rna_data))                     # features × cells
        X_h5, _, _ = fetch_matrix(src, nothing)                          # cells × features
        @test size(X_bp) == (462, 250)
        @test X_bp ≈ permutedims(X_h5) atol = 1e-4
        C_bp = Matrix(SparseArrays.sparse(rna_counts))
        C_h5, _, _ = fetch_matrix(src, nothing; layer = "counts")
        @test C_bp ≈ permutedims(C_h5) atol = 1e-6

        # as a modality inside MemoryAnnData, with obs/reductions from the h5ad
        rna, cells_ = SCUBA.bpcells_modality(Dict("data" => joinpath(dir, "RNA_data"), "counts" => joinpath(dir, "RNA_counts")))
        ab, _ = SCUBA.bpcells_modality(Dict("data" => joinpath(dir, "AB_data"), "counts" => joinpath(dir, "AB_counts")))
        bp = MemoryAnnData(; cells = cells_, obs = SCUBA.drop_cells(fetch_metadata(src)),
                           modalities = OrderedDict("RNA" => rna, "AB" => ab),
                           reductions = OrderedDict("UMAP" => Matrix(SCUBA.drop_cells(fetch_reduction(src)))),
                           default_layer = "data")
        test_source_conformance(bp; main_feature = "GAPDH", alt_modality = "AB", alt_feature = "CD117-AB",
                                reduction = "UMAP", categorical_var = "condensed_cell_type", numeric_var = "nCount_RNA")
        a = fetch_data(bp, ["RNA_GAPDH", "AB_CD117-AB", "GAPDH"]; layer = "data")
        b = fetch_data(src, ["X_GAPDH", "protein_CD117-AB", "GAPDH"])
        @test a.RNA_GAPDH ≈ b.X_GAPDH atol = 1e-4
        @test a[!, "AB_CD117-AB"] ≈ b[!, "protein_CD117-AB"] atol = 1e-4
        sub = all_cells(src)[[250, 3, 77]]
        @test fetch_data(bp, ["GAPDH"]; cells = sub).GAPDH ≈ fetch_data(src, ["GAPDH"]; cells = sub).GAPDH atol = 1e-4
        close(src)
    end
end

@testset "BPCells column construction" begin
    dir = fixture_path("AML_BPCells")
    if isdir(dir)
        m = open_bpcells(joinpath(dir, "RNA_data"))
        fidx = [5, 1, 300]; cidx = [10, 250, 1]
        @test reduce(hcat, SCUBA.pull_columns(m, fidx, cidx, Float64)) == SCUBA.pull(m, fidx, cidx)
        @test eltype(SCUBA.pull_columns(m, fidx, cidx, Float32)[1]) == Float32
    end
end
