# Build a Julia system image with SCUBA and its bridges compiled in, for fast cold starts from
# R (JuliaCall), Python (juliacall) or Julia itself.
#
#   julia --project=tools/sysimage tools/build_sysimage.jl [--out PATH] [--no-rcall] [--no-pythoncall]
#
# The default output is ~/.julia/scuba/SCUBA.<os>-<arch>.so (or .dylib / .dll). Use it with
#   R:      options(SCUBA.sysimage = "<path>")   before the first SCUBA call (or scuba_build_sysimage())
#   Python: PYTHON_JULIACALL_SYSIMAGE=<path>
#   Julia:  julia -J <path>
using Pkg
Pkg.activate(joinpath(@__DIR__, "sysimage"))
Pkg.develop(path = joinpath(@__DIR__, ".."); io = devnull)
# PackageCompiler needs every image package in the active project (Pkg.add is idempotent)
Pkg.add(["PackageCompiler", "RCall", "PythonCall", "DataFrames", "CategoricalArrays", "SparseArrays", "OrderedCollections", "Zarr"]; io = devnull)
Pkg.instantiate(io = devnull)
using PackageCompiler

args = copy(ARGS)
function flag!(name, default = nothing)
    i = findfirst(==(name), args)
    i === nothing && return default
    v = default === nothing ? args[i + 1] : true
    deleteat!(args, default === nothing ? (i:i + 1) : (i:i))
    return v
end
no_rcall = flag!("--no-rcall", false)
no_py = flag!("--no-pythoncall", false)
ext = Sys.isapple() ? ".dylib" : Sys.iswindows() ? ".dll" : ".so"
default_out = joinpath(homedir(), ".julia", "scuba", "SCUBA.$(Sys.KERNEL)-$(Sys.ARCH)$ext")
out = flag!("--out", default_out)
mkpath(dirname(out))

packages = [:SCUBA, :DataFrames, :CategoricalArrays, :Zarr]
no_rcall || push!(packages, :RCall)
no_py || push!(packages, :PythonCall)

workload = joinpath(@__DIR__, "sysimage_workload.jl")
@info "building sysimage" out packages
create_sysimage(packages; sysimage_path = out, precompile_execution_file = workload,
                include_transitive_dependencies = true)
println(out)
