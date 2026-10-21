# Build a Julia system image with SCUBA and its bridges compiled in, for fast cold starts from
# R (JuliaCall), Python (juliacall) or Julia itself.
#
#   julia tools/build_sysimage.jl [--out PATH] [--no-rcall] [--no-pythoncall] [--no-workload] [--cpu-target TARGET]
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
function switch!(name)              # boolean flag
    i = findfirst(==(name), args)
    i === nothing && return false
    deleteat!(args, i)
    return true
end
function option!(name, default)     # flag with a value
    i = findfirst(==(name), args)
    i === nothing && return default
    v = args[i + 1]
    deleteat!(args, i:i + 1)
    return v
end
no_rcall = switch!("--no-rcall")
no_py = switch!("--no-pythoncall")
no_workload = switch!("--no-workload")
cpu_target = option!("--cpu-target", PackageCompiler.default_app_cpu_target())
ext = Sys.isapple() ? ".dylib" : Sys.iswindows() ? ".dll" : ".so"
default_out = joinpath(homedir(), ".julia", "scuba", "SCUBA.$(Sys.KERNEL)-$(Sys.ARCH)$ext")
out = option!("--out", default_out)
isempty(args) || error("unknown arguments: $(join(args, ' '))")
mkpath(dirname(out))

packages = [:SCUBA, :DataFrames, :CategoricalArrays, :Zarr]
no_rcall || push!(packages, :RCall)
no_py || push!(packages, :PythonCall)

workload = no_workload ? String[] : joinpath(@__DIR__, "sysimage_workload.jl")
@info "building sysimage" out packages cpu_target
create_sysimage(packages; sysimage_path = out, precompile_execution_file = workload,
                include_transitive_dependencies = true, cpu_target)
println(out)
