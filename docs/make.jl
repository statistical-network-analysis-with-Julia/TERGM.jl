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
    repo = Documenter.Remotes.GitHub("Statistical-network-analysis-with-Julia", "TERGM.jl"),
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
    # STRICT. Undefined bindings, bad cross-references, duplicate docs and
    # malformed markdown are build ERRORS, so they cannot silently accumulate
    # again (a docs build that passes while warning is one that will rot).
    #
    # `checkdocs = :exports` is the one deliberate exclusion: every *exported*
    # name must be documented, but internal machinery (materialized/private
    # types, `Base`/`Graphs` method extensions, inner constructors) need not be
    # -- filler docstrings for names a user never types are worse than none.
    warnonly = false,
    checkdocs = :exports,
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
