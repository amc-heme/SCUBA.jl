# Interface conformance: invariants every source must satisfy. Called with a constructed source
# and the names of one main feature, one alt-modality feature and one reduction.
function test_source_conformance(src; main_feature, alt_modality, alt_feature, reduction, categorical_var, numeric_var)
    @testset "conformance: $(SCUBA.source_name(src))" begin
        n = length(all_cells(src))
        @test n > 0
        @test length(unique(all_cells(src))) == n
        mods = SCUBA.modality_keys(src)
        @test default_modality(src) == first(mods)
        @test alt_modality in mods
        main = default_modality(src)
        @test main_feature in features_in_modality(src)
        @test alt_feature in features_in_modality(src, alt_modality)
        @test default_layer(src) in SCUBA.layer_names(src, main)
        @test reduction in SCUBA.reduction_keys(src)
        @test default_reduction(src) == reduction
        @test categorical_var in meta_varnames(src)
        @test numeric_var in meta_varnames(src)

        keyed_main = string(main, "_", main_feature)
        keyed_alt = string(alt_modality, "_", alt_feature)
        red1 = SCUBA.reduction_dimnames(src, reduction, [1])[1]
        vars = [keyed_main, keyed_alt, red1, categorical_var, numeric_var]
        df = fetch_data(src, vars)
        @test names(df) == vcat(["cell"], vars)
        @test nrow(df) == n
        @test df.cell == all_cells(src)
        @test eltype(df[!, keyed_main]) == Float64
        @test eltype(df[!, keyed_alt]) == Float64
        @test eltype(df[!, red1]) == Float64
        @test SCUBA._is_categorical(df[!, categorical_var])
        @test nonmissingtype(eltype(df[!, numeric_var])) <: Real

        # unkeyed main feature equals keyed
        @test fetch_data(src, [main_feature])[!, main_feature] == df[!, keyed_main]
        # unkeyed alt feature is renamed with a warning and equals keyed
        (dfa, w) = with_collected_warnings(() -> fetch_data(src, [alt_feature]))
        @test names(dfa) == ["cell", keyed_alt]
        @test dfa[!, keyed_alt] == df[!, keyed_alt]
        @test any(x -> x.code == :alt_modality_renamed, w)
        # metadata + unkeyed main feature is silent
        (_, w2) = with_collected_warnings(() -> fetch_data(src, [main_feature, numeric_var]))
        @test isempty(w2)
        # cell subset preserves order and values
        sub = all_cells(src)[[3, 1, 2]]
        dfs = fetch_data(src, vars; cells = sub)
        @test dfs.cell == sub
        @test dfs[!, keyed_main] == df[[3, 1, 2], keyed_main]
        @test isequal(dfs[!, categorical_var], df[[3, 1, 2], categorical_var])
        # fetch_feature / fetch_matrix agree with fetch_data
        ff = fetch_feature(src, [main_feature, keyed_main])
        @test names(ff) == ["cell", main_feature]
        @test ff[!, main_feature] == df[!, keyed_main]
        # keyed names from another modality are fetched too, in the requested order, keeping their names
        ffm = fetch_feature(src, [keyed_alt, main_feature])
        @test names(ffm) == ["cell", keyed_alt, main_feature]
        @test ffm[!, keyed_alt] == df[!, keyed_alt] && ffm[!, main_feature] == df[!, keyed_main]
        @test_throws ArgumentError fetch_matrix(src, [keyed_alt, main_feature])
        # a bare alt-modality feature is not searched outside the main modality
        @test_throws SCUBA.VarsNotFound fetch_feature(src, [alt_feature])
        @test names(fetch_feature(src, [alt_feature]; modality = alt_modality)) == ["cell", alt_feature]
        X, c, f = fetch_matrix(src, [main_feature])
        @test vec(X) == df[!, keyed_main]
        Xs, _, _ = fetch_matrix(src, [main_feature]; sparse = true)
        @test Xs isa SparseMatrixCSC && Vector(vec(Xs)) == df[!, keyed_main]
        # sparse results with a requested element type agree with dense ones, also for reordered cells
        Xs32, _, _ = fetch_matrix(src, [main_feature]; sparse = true, eltype = Float32)
        @test Xs32 isa SparseMatrixCSC{Float32,Int} && size(Xs32) == (n, 1)
        @test Matrix(Xs32)[:, 1] ≈ Float32.(df[!, keyed_main])
        Xa32, _, _ = fetch_matrix(src, [alt_feature]; modality = alt_modality, sparse = true, eltype = Float32)
        @test Matrix(Xa32)[:, 1] ≈ Float32.(df[!, keyed_alt])
        @test Matrix(fetch_matrix(src, [main_feature]; cells = sub, sparse = true)[1]) == fetch_matrix(src, [main_feature]; cells = sub)[1]
        dsp = fetch_data(src, vars; sparse = true)
        @test names(dsp) == names(df)
        @test dsp[!, keyed_main] isa SparseVector && Vector(dsp[!, keyed_main]) == df[!, keyed_main]
        @test dsp[!, keyed_alt] isa SparseVector && Vector(dsp[!, keyed_alt]) == df[!, keyed_alt]
        @test isequal(dsp[!, categorical_var], df[!, categorical_var]) && dsp[!, red1] == df[!, red1]
        @test fetch_feature(src, [main_feature]; sparse = true)[!, main_feature] isa SparseVector
        @test_throws ArgumentError fetch_data(src, [keyed_main]; sparse = true, eltype = Int)
        # chunked iteration concatenates to the full result and keeps the requested cell order
        it = eachchunk(src, vars; cells_per_chunk = 2)
        @test length(it) == cld(n, 2)
        chunks = collect(it)
        @test all(nrow(c) <= 2 for c in chunks) && sum(nrow, chunks) == n
        cat = reduce(vcat, chunks)
        for v in names(df)
            @test isequal(cat[!, v], df[!, v])
        end
        catsub = reduce(vcat, collect(eachchunk(src, [keyed_main, numeric_var]; cells = sub, cells_per_chunk = 2)))
        @test catsub.cell == sub && catsub[!, keyed_main] == df[[3, 1, 2], keyed_main]
        @test isequal(catsub[!, numeric_var], df[[3, 1, 2], numeric_var])
        (_, wch) = with_collected_warnings(() -> collect(eachchunk(src, [alt_feature]; cells_per_chunk = 2)))
        @test count(x -> x.code == :alt_modality_renamed, wch) == 1
        @test reduce(vcat, collect(eachchunk(src, [keyed_main]; sparse = true, cells_per_chunk = 4)))[!, keyed_main] isa AbstractVector
        # metadata
        md = fetch_metadata(src, [categorical_var, numeric_var])
        @test names(md) == ["cell", categorical_var, numeric_var]
        @test isequal(md[!, categorical_var], df[!, categorical_var])
        @test isequal(fetch_metadata(src, numeric_var), df[!, numeric_var])
        full = fetch_metadata(src)
        @test names(full) == vcat(["cell"], meta_varnames(src))
        @test_throws MetadataNotFound fetch_metadata(src, [categorical_var, "__nope__"])
        # reductions
        fr = fetch_reduction(src)
        @test names(fr) == vcat(["cell"], SCUBA.reduction_dimnames(src, reduction, [1, 2]))
        @test fr[!, red1] == df[!, red1]
        @test_throws InvalidDims fetch_reduction(src; dims = [1, 10_000])
        @test_throws ReductionNotFound fetch_reduction(src, "__nope__")
        # cells / unique values
        uv = unique_values(src, categorical_var)
        @test length(uv) == n_unique(src, categorical_var)
        lvl = string(first(uv))
        cs = fetch_cells(src, categorical_var, [lvl])
        @test !isempty(cs) && issubset(cs, all_cells(src))
        @test cs == df.cell[string.(df[!, categorical_var]) .== lvl]
        @test_throws NonCategoricalVar fetch_cells(src, numeric_var, [1])
        @test_throws NonCategoricalVar unique_values(src, numeric_var)
        # errors
        @test_throws VarsNotFound fetch_data(src, ["__definitely_not_a_var__"])
        @test_throws CellsNotFound fetch_data(src, [keyed_main]; cells = ["__nope__"])
        # snapshot round trip
        snap = snapshot(src)
        @test isequal(fetch_data(snap, vars), df)
        @test all_keys(snap) == all_keys(src)
    end
end
