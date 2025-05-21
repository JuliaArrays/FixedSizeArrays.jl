using Documenter, FixedSizeArrays

makedocs(
    modules = [FixedSizeArrays],
    sitename = "FixedSizeArrays",
    pages    = [
        "Introduction" => "index.md",
        "Installation" => "installation.md",
        "Usage" => "usage.md",
        "Reference" => "reference.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaArrays/FixedSizeArrays.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
    push_preview = true,
)
