# Compare SCUBA.jl (h5ad backend) against golden outputs recorded from SCUBA R 1.2.2 on the
# SingleCellExperiment form of the same data (fixtures/build/goldens, produced by make_goldens.R).
# Column names differ by form (RNA_ vs X_, UMAP_ vs X_umap_), so columns are compared by position.
using DelimitedFiles

const GOLDEN_DIR = let d = get(ENV, "SCUBA_GOLDEN_DIR", "")
    isempty(d) ? normpath(joinpath(@__DIR__, "..", "fixtures", "build", "goldens")) : d
end

function read_golden_csv(path)
    raw, header = readdlm(path, ','; header = true, quotes = true)
    cols = String.(vec(header))
    data = Dict{String,Vector{Any}}()
    for (j, c) in enumerate(cols)
        data[c] = Any[raw[i, j] for i in 1:size(raw, 1)]
    end
    return cols, data
end

_num(x) = x isa Number ? Float64(x) : (isempty(String(x)) ? NaN : parse(Float64, String(x)))
_str(x) = x isa Number ? string(x) : String(x)
_jlnum(x) = ismissing(x) ? NaN : Float64(x)
_jlstr(x) = ismissing(x) ? "" : string(x)
_approx_or_both_nan(a, b; atol) = (isnan(a) && isnan(b)) || isapprox(a, b; atol)

function compare_numeric(jl::AbstractVector, gold::Vector{Any}; atol = 1e-6)
    @test length(jl) == length(gold)
    @test all(_approx_or_both_nan.(_jlnum.(jl), _num.(gold); atol))
end
function compare_categorical(jl::AbstractVector, gold::Vector{Any})
    @test length(jl) == length(gold)
    @test _jlstr.(jl) == _str.(gold)
end

@testset "goldens vs SCUBA R 1.2.2 (SCE form)" begin
    h5 = fixture_path("AML_h5ad.h5ad")
    if !(isdir(GOLDEN_DIR) && isfile(h5))
        @warn "goldens or h5ad fixture not found; skipping golden tests" GOLDEN_DIR
    else
        src = open_h5ad(h5)
        gold(id) = read_golden_csv(joinpath(GOLDEN_DIR, "$(id)__sce.csv"))
        goldtxt(id) = readlines(joinpath(GOLDEN_DIR, "$(id)__sce.txt"))

        @testset "basic_mixed" begin
            cols, g = gold("basic_mixed")
            df = fetch_data(src, ["X_GAPDH", "protein_CD117-AB", "X_umap_1", "X_umap_2", "nCount_RNA", "condensed_cell_type", "Batch"])
            @test df.cell == String.(g["cell"])
            @test length(names(df)) == length(cols)
            for j in 2:5
                compare_numeric(df[!, j], g[cols[j]])
            end
            compare_numeric(df.nCount_RNA, g["nCount_RNA"])
            compare_categorical(df.condensed_cell_type, g["condensed_cell_type"])
            compare_categorical(df.Batch, g["Batch"])
        end

        @testset "unkeyed main + metadata (silent)" begin
            cols, g = gold("unkeyed_main_and_metadata")
            (df, w) = with_collected_warnings(() -> fetch_data(src, ["GAPDH", "MEIS1", "nCount_RNA"]))
            @test isempty(w)
            @test names(df)[2:end] == cols[2:end]
            for c in cols[2:end]
                compare_numeric(df[!, c], g[c])
            end
        end

        @testset "unkeyed alt feature renamed with warning" begin
            cols, g = gold("unkeyed_alt_feature")
            @test cols[2] == "AB_CD11a-AB"
            (df, w) = with_collected_warnings(() -> fetch_data(src, ["CD11a-AB"]))
            @test names(df)[2] == "protein_CD11a-AB"
            @test any(x -> x.code == :alt_modality_renamed, w)
            compare_numeric(df[!, 2], g[cols[2]])
        end

        @testset "partial missing: warning, found columns returned" begin
            cols, g = gold("partial_missing")
            @test cols[2:end] == ["RNA_GAPDH"]
            (df, w) = with_collected_warnings(() -> fetch_data(src, ["X_GAPDH", "nonsense", "X_MESS1", "protein_CD1110-AB", "X_umap_11"]))
            @test names(df)[2:end] == ["X_GAPDH"]
            @test any(x -> x.code == :some_vars_not_found && x.vars == ["nonsense", "X_MESS1", "protein_CD1110-AB", "X_umap_11"], w)
            compare_numeric(df.X_GAPDH, g["RNA_GAPDH"])
        end

        @testset "all missing: error" begin
            @test !isfile(joinpath(GOLDEN_DIR, "all_missing__sce.csv"))
            @test_throws VarsNotFound fetch_data(src, ["protein_CD900"])
        end

        @testset "duplicates" begin
            cols, g = gold("duplicates_literal")
            (df, w) = with_collected_warnings(() -> fetch_data(src, ["X_GAPDH", "X_GAPDH", "X_MEIS1", "X_MEIS1"]))
            @test length(names(df)) == length(cols)
            @test any(x -> x.code == :duplicate_var, w)
            cols2, g2 = gold("duplicates_key_equivalent")
            @test cols2[2:end] == ["RNA_GAPDH", "RNA_MEIS1"]
            (df2, w2) = with_collected_warnings(() -> fetch_data(src, ["X_GAPDH", "GAPDH", "X_MEIS1", "MEIS1"]))
            @test names(df2)[2:end] == ["X_GAPDH", "X_MEIS1"]
            compare_numeric(df2.X_GAPDH, g2["RNA_GAPDH"])
        end

        @testset "counts layer" begin
            cols, g = gold("counts_layer")
            df = fetch_data(src, ["X_GAPDH", "X_UNG"]; layer = "counts")
            compare_numeric(df.X_GAPDH, g["RNA_GAPDH"])
            compare_numeric(df.X_UNG, g["RNA_UNG"])
        end

        @testset "cells subset" begin
            # legacy divergence: SCUBA R 1.2.2 mishandles a single reduction dimension on SCE objects
            # (the column comes back named `data`), so this case is checked against the Seurat golden.
            cols, g = read_golden_csv(joinpath(GOLDEN_DIR, "cells_subset__seurat.csv"))
            @test cols == ["cell", "rna_GAPDH", "UMAP_1", "Batch"]
            sub = reverse(all_cells(src)[1:20])
            df = fetch_data(src, ["X_GAPDH", "X_umap_1", "Batch"]; cells = sub)
            @test df.cell == String.(g["cell"])
            compare_numeric(df.X_GAPDH, g["rna_GAPDH"])
            compare_numeric(df.X_umap_1, g["UMAP_1"])
            compare_categorical(df.Batch, g["Batch"])
        end

        @testset "fetch_reduction default" begin
            cols, g = gold("reduction_default")
            df = fetch_reduction(src)
            compare_numeric(df[!, 2], g[cols[2]])
            compare_numeric(df[!, 3], g[cols[3]])
        end

        @testset "fetch_metadata" begin
            cols, g = gold("metadata_vars")
            df = fetch_metadata(src, ["nCount_RNA", "ct", "prop"])
            compare_numeric(df.nCount_RNA, g["nCount_RNA"])
            compare_categorical(df.ct, g["ct"])
            compare_numeric(df.prop, g["prop"])
            @test_throws MetadataNotFound fetch_metadata(src, ["nCount_RNA", "nonsense"])
            colsf, gf = gold("metadata_full")
            full = fetch_metadata(src)
            # the SCE fixture carries two extra columns (orig.ident, ident) that the h5ad lacks
            shared = intersect(names(full), colsf)
            @test length(shared) >= 28
            for c in shared
                c == "cell" && continue
                if SCUBA._is_categorical(full[!, c])
                    compare_categorical(full[!, c], gf[c])
                else
                    compare_numeric(full[!, c], gf[c])
                end
            end
        end

        @testset "fetch_cells / unique_values / meta_varnames" begin
            @test fetch_cells(src, "Batch", ["PBMC_200AB"]) == goldtxt("cells_by_batch")
            @test_throws NonCategoricalVar fetch_cells(src, "nCount_RNA", [1])
            @test string.(unique_values(src, "condensed_cell_type")) == goldtxt("unique_cell_types")
            mv = goldtxt("meta_varnames")
            @test issubset(meta_varnames(src), mv)
            @test filter(x -> x in meta_varnames(src), mv) == meta_varnames(src)   # same order
        end
        close(src)
    end
end
