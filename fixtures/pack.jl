# Pack fixtures/data into a tarball for the lazy `scuba_fixtures` artifact and print the
# Artifacts.toml entry. Upload the tarball to a GitHub release, then paste the entry (with the
# release URL) into Artifacts.toml.
#
#   julia --project=. fixtures/pack.jl [url-prefix]
using Pkg.Artifacts, SHA, Tar, CodecZlib
data = joinpath(@__DIR__, "data")
isdir(data) || error("fixtures/data not found")
version = "v1"
out = joinpath(@__DIR__, "build", "scuba-fixtures-$version.tar.gz")
mkpath(dirname(out))
open(out, "w") do io
    gz = GzipCompressorStream(io)
    Tar.create(data, gz)
    close(gz)
end
tree = Pkg.Artifacts.create_artifact() do dir
    cp(data, dir; force = true)
end
sha = bytes2hex(open(sha256, out))
prefix = length(ARGS) >= 1 ? ARGS[1] : "https://github.com/amc-heme/SCUBA.jl/releases/download/fixtures-$version"
println("""
[scuba_fixtures]
git-tree-sha1 = "$tree"
lazy = true

    [[scuba_fixtures.download]]
    url = "$prefix/$(basename(out))"
    sha256 = "$sha"
""")
println("tarball: $out  ($(round(filesize(out) / 2^20; digits = 1)) MiB)")
