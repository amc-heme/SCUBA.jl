# Execution workload recorded while building the sysimage (PackageCompiler
# `precompile_execution_file`). Exercises the public API on in-memory and on-disk sources; the
# RCall / PythonCall bridges are loaded but not exercised (they need a live host runtime).
using SCUBA, DataFrames, CategoricalArrays, SparseArrays, OrderedCollections
SCUBA._precompile_workload()
try
    path = SCUBA.fixture_path("AML_h5ad.h5ad")
    if isfile(path)
        src = open_h5ad(path)
        fetch_data(src, ["X_GAPDH", "protein_CD117-AB", "X_umap_1", "Batch", "nCount_RNA"])
        fetch_data(src, ["GAPDH"]; cells = all_cells(src)[1:10], layer = "counts", eltype = Float32)
        fetch_metadata(src); fetch_reduction(src); fetch_cells(src, "Batch", ["PBMC_200AB"]); unique_values(src, "condensed_cell_type")
        snapshot(src); close(src)
    end
    bp = SCUBA.fixture_path("AML_BPCells/RNA_data")
    if isdir(bp)
        m = open_bpcells(bp); SCUBA.pull(m, [1, 2], [1, 2, 3]); SCUBA.pull_columns(m, [1], [1, 2], Float64)
    end
catch err
    @warn "sysimage workload: fixture part skipped" err
end
