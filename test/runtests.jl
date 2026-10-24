using Test
using SCUBA
using DataFrames, CategoricalArrays, SparseArrays, OrderedCollections

# Pass "bridges" (Pkg.test(test_args=["bridges"])) to also run the R / Python / Zarr extension tests.
# They need: Zarr.jl; RCall.jl with an R that has SeuratObject, SingleCellExperiment, HDF5Array (and
# BPCells for the BPCells case); PythonCall.jl with JULIA_PYTHONCALL_EXE pointing at a Python with anndata.
const RUN_BRIDGES = "bridges" in ARGS
const ONLY = filter(a -> startswith(a, "only="), ARGS)
const BRIDGE_SET = isempty(ONLY) ? ["zarr", "rcall", "pythoncall"] : split(replace(ONLY[1], "only=" => ""), ",")

include("synthetic.jl")
include("conformance.jl")
include("test_resolve.jl")
include("test_h5ad.jl")
include("test_h5mu.jl")
include("test_memory.jl")
include("test_overlay.jl")
include("test_bpcells.jl")
include("test_aggregate.jl")
include("test_golden.jl")

if RUN_BRIDGES
    # RCall / PythonCall are not part of the test target (installing RCall without R falls back to
    # Conda, which is slow and breaks on nightly). Add the requested bridges to the sandboxed test
    # environment on demand; SCUBA's [compat] bounds still apply.
    using Pkg
    let want = String[]
        "rcall" in BRIDGE_SET && push!(want, "RCall")
        "pythoncall" in BRIDGE_SET && push!(want, "PythonCall")
        isempty(want) || Pkg.add(want)
    end
    "zarr" in BRIDGE_SET && include("bridges/test_zarr.jl")
    "pythoncall" in BRIDGE_SET && include("bridges/test_pythoncall.jl")
    "rcall" in BRIDGE_SET && include("bridges/test_rcall.jl")
else
    @info "Bridge tests skipped (pass test_args=[\"bridges\"] to run them)"
end
