using Documenter, SCUBA

makedocs(;
    modules = [SCUBA],
    sitename = "SCUBA.jl",
    authors = "SCUBA Authors",
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true",
                             canonical = "https://amc-heme.github.io/SCUBA.jl"),
    pages = ["Home" => "index.md",
             "Fetching data" => "fetching.md",
             "Backends" => "backends.md",
             "Compatibility with SCUBA (R)" => "compat.md",
             "API reference" => "api.md"],
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(; repo = "github.com/amc-heme/SCUBA.jl", devbranch = "main", push_preview = true)
