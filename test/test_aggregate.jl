# aggregate(): group sums / means of a layer without materialising the cells × features block.

# A lazy layer for the in-memory backend (values doubled on read), defined at top level because
# struct definitions are not allowed inside a testset block.
struct AggregateDoublingLazy <: SCUBA.LazyMatrix
    A::Matrix{Float64}
end
Base.size(m::AggregateDoublingLazy) = size(m.A)
SCUBA.pull(m::AggregateDoublingLazy, fidx, cidx) = 2 .* m.A[cidx, fidx]

# cells × groups indicator matrix for a code vector (0 = cell left out)
function aggregate_indicator(codes, ngroups)
    M = zeros(Float64, length(codes), ngroups)
    for (i, g) in enumerate(codes)
        g == 0 || (M[i, g] = 1.0)
    end
    return M
end
# group codes (1-based, sorted level order) of a categorical / string column
function aggregate_codes(col)
    s = string.(col)
    lv = sort(unique(s))
    return Int32[findfirst(==(x), lv) for x in s], length(lv)
end
# reference: indicator' * dense(cells × features) read through fetch_matrix
function aggregate_reference(src, codes, ngroups; kw...)
    X = Matrix(fetch_matrix(src, nothing; sparse = true, kw...)[1])
    return aggregate_indicator(codes, ngroups)' * X
end

@testset "aggregate on the synthetic source" begin
    src = synthetic_source()
    cells = all_cells(src)
    counts = Matrix(fetch_matrix(src, nothing; layer = "counts")[1])      # 6 cells × 5 features
    codes = Int32[1, 1, 2, 2, 3, 0]
    V, n, fs, cs = aggregate(src, codes, 3; layer = "counts")
    @test V isa Matrix{Float64} && size(V) == (3, 5)
    @test V == aggregate_indicator(codes, 3)' * counts
    @test n == [2, 2, 1]
    @test fs == features_in_modality(src) && cs == cells
    # mean = sum / contributing cells
    Vm, nm, = aggregate(src, codes, 3; layer = "counts", fun = :mean)
    @test Vm ≈ V ./ [2, 2, 1] && nm == n
    # an empty group stays zero (and is not divided by zero for the mean)
    V0, n0, = aggregate(src, Int32[1, 1, 0, 0, 0, 0], 2; layer = "counts", fun = :mean)
    @test n0 == [2, 0] && all(iszero, V0[2, :]) && V0[1, :] == (counts[1, :] .+ counts[2, :]) ./ 2
    # plain Int codes, a reordered cell subset
    sub = cells[[5, 2, 6]]
    Vs, ns, _, css = aggregate(src, [1, 2, 1], 2; cells = sub, layer = "counts")
    @test css == sub && ns == [2, 1]
    @test Vs == permutedims([counts[5, :] .+ counts[6, :] counts[2, :]])
    # feature subset keeps the requested order; keyed names resolve
    Vf, _, fsf = aggregate(src, codes, 3; layer = "counts", features = ["g2", "RNA_g1"])
    @test fsf == ["g2", "RNA_g1"] && Vf == V[:, [2, 1]]
    # default layer when none is given
    @test aggregate(src, codes, 3)[1] == aggregate_indicator(codes, 3)' * Matrix(fetch_matrix(src, nothing)[1])
    # a dense alternative modality
    ab = Matrix(fetch_matrix(src, nothing; modality = "AB")[1])
    Va, na, fsa = aggregate(src, codes, 3; modality = "AB")
    @test Va == aggregate_indicator(codes, 3)' * ab && na == [2, 2, 1] && fsa == features_in_modality(src, "AB")
    # chunking never changes the result
    for k in (1, 2, 4, 10)
        @test aggregate(src, codes, 3; layer = "counts", cells_per_chunk = k)[1] == V
    end
    # the sink itself
    s = SCUBA.GroupSumSink(Float64, codes, 3, 5)
    SCUBA.expression!(s, src, "RNA", "counts", collect(1:5), collect(1:6))
    @test s.acc == V
    # lazy layers go through pull!
    A = Float64[i * 10 + j for i in 1:4, j in 1:3]
    lazy = MemoryAnnData(; cells = ["a", "b", "c", "d"],
                         modalities = OrderedDict("M" => Modality(["f1", "f2", "f3"], OrderedDict("X" => AggregateDoublingLazy(A)))))
    @test aggregate(lazy, [1, 2, 1, 2], 2)[1] == aggregate_indicator([1, 2, 1, 2], 2)' * (2 .* A)
    # the overlay forwards
    ov = MetadataOverlay(src; obs = DataFrame())
    @test aggregate(ov, codes, 3; layer = "counts")[1] == V

    @testset "argument errors" begin
        @test_throws ArgumentError aggregate(src, Int32[1, 2], 2; layer = "counts")             # wrong length
        @test_throws ArgumentError aggregate(src, Int32[1, 1, 2, 2, 3, 4], 3; layer = "counts")  # code > ngroups
        @test_throws ArgumentError aggregate(src, Int32[-1, 1, 2, 2, 3, 0], 3; layer = "counts") # negative code
        @test_throws ArgumentError aggregate(src, codes, 3; layer = "counts", fun = :max)
        @test_throws ArgumentError aggregate(src, codes, -1; layer = "counts")
        @test_throws ArgumentError aggregate(src, codes, 3; features = ["g1", "AB_p1"])          # two modalities
        @test_throws SCUBA.LayerNotFound aggregate(src, codes, 3; layer = "nope")
        @test_throws SCUBA.CellsNotFound aggregate(src, [1], 1; cells = ["zz"])
    end
end

@testset "aggregate on h5ad (csc, csr, dense) and the protein data frame" begin
    paths = [fixture_path(f) for f in ("AML_h5ad.h5ad", "AML_h5ad_csr.h5ad", "AML_h5ad_dense.h5ad")]
    if !all(isfile, paths)
        @warn "h5ad fixtures not found; skipping aggregate h5ad tests"
    else
        base = open_h5ad(paths[1])
        codes, ng = aggregate_codes(fetch_metadata(base, "condensed_cell_type"))
        cells = all_cells(base)
        sub = cells[[200, 3, 77, 4, 250, 1, 150]]
        for p in paths
            s = open_h5ad(p)
            for layer in SCUBA.layer_names(s, "X")
                want = aggregate_reference(s, codes, ng; layer = layer)
                V, n, fs, = aggregate(s, codes, ng; layer = layer)
                @test V ≈ want && n == vec(sum(aggregate_indicator(codes, ng); dims = 1)) && fs == features_in_modality(s)
                # bounded reads and explicit chunks give the same sums
                old = SCUBA.MAX_READ_ELEMENTS[]
                SCUBA.MAX_READ_ELEMENTS[] = 100
                try
                    @test aggregate(s, codes, ng; layer = layer)[1] ≈ want
                    @test aggregate(s, codes, ng; layer = layer, cells_per_chunk = 7)[1] ≈ want
                finally
                    SCUBA.MAX_READ_ELEMENTS[] = old
                end
                # a reordered cell subset
                Xs = Matrix(fetch_matrix(s, nothing; layer = layer, cells = sub, sparse = true)[1])
                Vs, ns, = aggregate(s, [1, 2, 1, 2, 3, 3, 1], 3; layer = layer, cells = sub)
                @test Vs ≈ aggregate_indicator([1, 2, 1, 2, 3, 3, 1], 3)' * Xs && ns == [3, 2, 2]
            end
            # the data-frame-encoded protein modality
            Vp, np_, fsp = aggregate(s, codes, ng; modality = "protein")
            @test Vp ≈ aggregate_reference(s, codes, ng; modality = "protein") && fsp == features_in_modality(s, "protein")
            @test np_ == vec(sum(aggregate_indicator(codes, ng); dims = 1))
            close(s)
        end
        # ScubaR-style overlay over a file source
        ov = MetadataOverlay(base; obs = DataFrame())
        @test aggregate(ov, codes, ng; layer = "counts")[1] ≈ aggregate_reference(base, codes, ng; layer = "counts")
        close(base)
    end
end

@testset "aggregate on a BPCells directory" begin
    dir = fixture_path("AML_BPCells")
    h5 = fixture_path("AML_h5ad.h5ad")
    if !(isdir(dir) && isfile(h5))
        @warn "BPCells fixtures not found; skipping aggregate BPCells tests"
    else
        base = open_h5ad(h5)
        codes, ng = aggregate_codes(fetch_metadata(base, "condensed_cell_type"))
        want = aggregate_reference(base, codes, ng; layer = "counts")
        fsh = features_in_modality(base)
        # a BPCells directory is a layer; wrap it as the modality of an in-memory container
        rna, cells_ = SCUBA.bpcells_modality(Dict("counts" => joinpath(dir, "RNA_counts")))
        m = MemoryAnnData(; cells = cells_, modalities = OrderedDict("RNA" => rna), default_layer = "counts")
        V, n, fsb, = aggregate(m, codes, ng)
        @test n == vec(sum(aggregate_indicator(codes, ng); dims = 1))
        @test V[:, [findfirst(==(f), fsb) for f in fsh]] ≈ want
        @test aggregate(m, codes, ng; cells_per_chunk = 13)[1] == V
        sub = all_cells(m)[[10, 1, 99, 250]]
        Vs, = aggregate(m, [1, 1, 2, 2], 2; cells = sub)
        @test Vs ≈ aggregate_indicator([1, 1, 2, 2], 2)' * Matrix(fetch_matrix(m, nothing; cells = sub, sparse = true)[1])
        # the native reader streams: no chunking hook
        @test SCUBA.aggregate_chunk_cells(m, "RNA", "counts", 462) === nothing
        close(base)
    end
end

@testset "aggregate on MuData (full and partial membership)" begin
    h5 = fixture_path("AML_h5ad.h5ad")
    full = fixture_path("AML.h5mu")
    partial = fixture_path("AML_partial.h5mu")
    if !(isfile(h5) && isfile(full))
        @warn "MuData fixtures not found; skipping aggregate MuData tests"
    else
        base = open_h5ad(h5)
        codes, ng = aggregate_codes(fetch_metadata(base, "condensed_cell_type"))
        mu = open_h5mu(full)
        @test all_cells(mu) == all_cells(base)
        Vr, nr, = aggregate(mu, codes, ng; modality = "rna", layer = "counts")
        @test Vr ≈ aggregate_reference(base, codes, ng; layer = "counts") && nr == vec(sum(aggregate_indicator(codes, ng); dims = 1))
        Vp, = aggregate(mu, codes, ng; modality = "prot")
        @test Vp ≈ aggregate_reference(base, codes, ng; modality = "protein")
        close(mu)
        close(base)
    end
    # partial membership (fixture built by fixtures/make_mudata_partial.py; not in the artifact)
    if !isfile(partial)
        @warn "AML_partial.h5mu not found; skipping the partial-membership aggregate tests"
    else
        mp = open_h5mu(partial)
        codes_p, ng_p = aggregate_codes(fetch_metadata(mp, "condensed_cell_type"))
        absent = findall(==(0), mp.modalities["prot"].obsmap)
        @test length(absent) == 50
        Vq, nq, = aggregate(mp, codes_p, ng_p; modality = "prot")
        # absent cells contribute nothing and are not counted
        @test sum(nq) == 250 - 50
        @test nq == vec(sum(aggregate_indicator(codes_p, ng_p)[setdiff(1:250, absent), :]; dims = 1))
        X = Matrix(fetch_matrix(mp, nothing; modality = "prot")[1])      # NaN rows for absent cells
        X[isnan.(X)] .= 0.0
        @test Vq ≈ aggregate_indicator(codes_p, ng_p)' * X
        @test aggregate(mp, codes_p, ng_p; modality = "prot", fun = :mean)[1] ≈ Vq ./ max.(nq, 1)
        # a cell subset mixing absent and present cells
        sub = all_cells(mp)[[absent[1], 2, absent[2], 7]]
        Vsub, nsub, = aggregate(mp, [1, 1, 2, 2], 2; modality = "prot", cells = sub)
        @test nsub == [1, 1]
        Xsub = Matrix(fetch_matrix(mp, nothing; modality = "prot", cells = sub)[1]); Xsub[isnan.(Xsub)] .= 0.0
        @test Vsub ≈ aggregate_indicator([1, 1, 2, 2], 2)' * Xsub
        close(mp)
    end
end
