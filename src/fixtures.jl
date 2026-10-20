# Location of the shared test/demo fixtures (250-cell AML CITE-seq dataset in several formats).

const FIXTURES_ENV = "SCUBA_FIXTURES_DIR"

"""
    fixtures_dir() -> String

Directory holding the AML fixtures. Resolution order: the `SCUBA_FIXTURES_DIR` environment
variable, a `fixtures/data` directory in a development checkout, then the lazy
`scuba_fixtures` artifact (downloaded on first use).
"""
function fixtures_dir()
    env = get(ENV, FIXTURES_ENV, "")
    isempty(env) || isdir(env) || @warn "$FIXTURES_ENV is set but is not a directory; ignoring" env
    isempty(env) || !isdir(env) || return env
    dev = normpath(joinpath(@__DIR__, "..", "fixtures", "data"))
    isdir(dev) && isfile(joinpath(dev, "AML_h5ad.h5ad")) && return dev
    toml = joinpath(@__DIR__, "..", "Artifacts.toml")
    if isfile(toml) && artifact_hash("scuba_fixtures", toml) !== nothing
        LazyArtifacts.ensure_artifact_installed("scuba_fixtures", toml)
        return artifact_path(artifact_hash("scuba_fixtures", toml))
    end
    error("SCUBA fixtures not found. Set $FIXTURES_ENV to a directory containing AML_h5ad.h5ad etc.")
end

"""`fixture_path(name)` — absolute path of a fixture file or directory."""
fixture_path(name::AbstractString) = joinpath(fixtures_dir(), name)
