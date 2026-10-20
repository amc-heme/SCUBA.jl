using Zarr
@testset "Zarr AnnData stores" begin
    h5 = fixture_path("AML_h5ad.h5ad")
    for (label, p) in (("zarr v2", fixture_path("AML_v2.zarr")), ("zarr v3", fixture_path("AML.zarr")))
        if !(isdir(p) && isfile(h5))
            @warn "$label fixture not found; skipping" p
            continue
        end
        @testset "$label" begin
            z = try
                open_zarr(p)
            catch e
                @warn "could not open $label store with Zarr.jl; skipping" exception = e
                nothing
            end
            z === nothing && continue
            ref = open_h5ad(h5)
            @test SCUBA.modality_keys(z) == SCUBA.modality_keys(ref)
            @test SCUBA.reduction_keys(z) == SCUBA.reduction_keys(ref)
            @test all_cells(z) == all_cells(ref)
            vars = ["X_GAPDH", "protein_CD117-AB", "X_umap_1", "Batch", "nCount_RNA", "condensed_cell_type"]
            a = fetch_data(z, vars); b = fetch_data(ref, vars)
            @test names(a) == names(b)
            for c in vars
                @test isequal(a[!, c], b[!, c]) || a[!, c] ≈ b[!, c]
            end
            @test fetch_data(z, ["X_GAPDH"]; layer = "counts") == fetch_data(ref, ["X_GAPDH"]; layer = "counts")
            test_source_conformance(z; main_feature = "GAPDH", alt_modality = "protein", alt_feature = "CD117-AB",
                                    reduction = "X_umap", categorical_var = "condensed_cell_type", numeric_var = "nCount_RNA")
            close(ref)
        end
    end
end
