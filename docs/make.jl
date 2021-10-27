using Documenter
using NoiseCC


makedocs(
    modules = [NoiseCC],
    sitename = "NoiseCC",
    authors = "Fu Yin",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    pages = Any[
        "Home" => "index.md",
        "Framwork" => "Framwork.md",
        "FFT" => "FFT.md",
        "CC" => "CC.md",
        "FFT_CC" => "FFT_CC.md",
        "STACK" => "STACK.md",
        "Plot" => "plot.md",
        "References" => "references.md",
        ],
)

# deploydocs(
#     repo = "github.com/OUCyf/NoiseCC.jl.git",
#     target = "build",
#     deps   = nothing,
#     make   = nothing,
# )