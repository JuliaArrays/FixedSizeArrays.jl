using Documenter: Documenter, makedocs, deploydocs
using FixedSizeArrays
using Changelog: Changelog

const GITHUB_REPO = "JuliaArrays/FixedSizeArrays.jl"

# Generate a Documenter-friendly changelog from CHANGELOG.md
Changelog.generate(
    Changelog.Documenter(),
    joinpath(@__DIR__, "..", "CHANGELOG.md"),
    joinpath(@__DIR__, "src", "release-notes.md");
    repo = GITHUB_REPO,
    branch = "main",
)

makedocs(
    modules = [FixedSizeArrays],
    sitename = "FixedSizeArrays",
    pages    = [
        "Introduction" => "index.md",
        "Installation" => "installation.md",
        "Usage" => "usage.md",
        "Reference" => "reference.md",
        "Release Notes" => "release-notes.md",
    ],
    format = Documenter.HTML(
        ;
        prettyurls = haskey(ENV, "CI"),
    ),
)

deploydocs(
    repo = "github.com/$(GITHUB_REPO).git",
    target = "build",
    deps = nothing,
    make = nothing,
    push_preview = true,
)
