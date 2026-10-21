# Build a Julia system image with SCUBA (and optionally its bridges) compiled in, for fast cold
# starts from R (JuliaCall), Python (juliacall) or Julia itself. Needs a C compiler on the PATH
# (gcc/clang; `apt install build-essential` on Debian/Ubuntu) and network access for Pkg. On
# Julia 1.13 RCall cannot be included (LLVM error), so use --no-rcall for the R engine's image.
#
#   julia tools/build_sysimage.jl [--out PATH] [--no-rcall] [--no-pythoncall] [--zarr] [--no-workload] [--cpu-target TARGET]
#
# The default output is ~/.julia/scuba/SCUBA.<os>-<arch>.so (or .dylib / .dll). Use it with
#   R:      options(SCUBA.sysimage = "<path>")   before the first SCUBA call (or scuba_build_sysimage())
#   Python: PYTHON_JULIACALL_SYSIMAGE=<path>
#   Julia:  julia -J <path>
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
with_zarr = switch!("--zarr")
# "native" compiles one copy of the code for this machine's CPU. PackageCompiler's portable default
# compiles several CPU-target clones, which roughly triples compile time and memory (a build was
# OOM-killed on a server with the default); pass --cpu-target with PackageCompiler's default to
# get a portable image: "generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1)".
cpu_target = option!("--cpu-target", "native")

using Pkg
mkpath(joinpath(@__DIR__, "sysimage"))     # build project, created on demand (not tracked in git)
Pkg.activate(joinpath(@__DIR__, "sysimage"))
Pkg.develop(path = joinpath(@__DIR__, ".."); io = devnull)
# PackageCompiler needs every image package in the active project (Pkg.add is idempotent). Only
# the bridges that go into the image are added: RCall's build step needs an R installation and
# PythonCall drags in CondaPkg, so neither should be required for an image that excludes it.
wanted = ["PackageCompiler", "DataFrames", "CategoricalArrays", "SparseArrays", "OrderedCollections"]
no_rcall || push!(wanted, "RCall")
no_py || push!(wanted, "PythonCall")
with_zarr && push!(wanted, "Zarr")
Pkg.add(wanted; io = devnull)
for bridge in ("RCall", "PythonCall", "Zarr")   # left over from an earlier build with other flags
    bridge in wanted || !haskey(Pkg.project().dependencies, bridge) || Pkg.rm(bridge; io = devnull)
end
Pkg.instantiate(io = devnull)
using PackageCompiler

ext = Sys.isapple() ? ".dylib" : Sys.iswindows() ? ".dll" : ".so"
default_out = joinpath(homedir(), ".julia", "scuba", "SCUBA.$(Sys.KERNEL)-$(Sys.ARCH)$ext")
out = option!("--out", default_out)
isempty(args) || error("unknown arguments: $(join(args, ' '))")
mkpath(dirname(out))

packages = [:SCUBA, :DataFrames, :CategoricalArrays]
no_rcall || push!(packages, :RCall)
no_py || push!(packages, :PythonCall)
with_zarr && push!(packages, :Zarr)

workload = no_workload ? String[] : joinpath(@__DIR__, "sysimage_workload.jl")
@info "building sysimage" out packages cpu_target
create_sysimage(packages; sysimage_path = out, precompile_execution_file = workload,
                include_transitive_dependencies = true, cpu_target)
println(out)
