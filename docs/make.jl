using Documenter
using TERGM

DocMeta.setdocmeta!(TERGM, :DocTestSetup, :(using TERGM); recursive=true)

makedocs(
    sitename = "TERGM.jl",
    modules = [TERGM],
    authors = "Statistical Network Analysis with Julia",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://Statistical-network-analysis-with-Julia.github.io/TERGM.jl",
        edit_link = "main",
    ),
    repo = "https://github.com/Statistical-network-analysis-with-Julia/TERGM.jl/blob/{commit}{path}#{line}",
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "User Guide" => [
            "Temporal Terms" => "guide/terms.md",
            "STERGM Estimation" => "guide/estimation.md",
            "Simulation" => "guide/simulation.md",
        ],
        "API Reference" => [
            "Types" => "api/types.md",
            "Terms" => "api/terms.md",
            "Estimation" => "api/estimation.md",
        ],
    ],
    warnonly = [:missing_docs, :docs_block],
)

deploydocs(
    repo = "github.com/Statistical-network-analysis-with-Julia/TERGM.jl.git",
    devbranch = "main",
    versions = [
        "stable" => "dev",
        "dev" => "dev",
    ],
    push_preview = true,
)
