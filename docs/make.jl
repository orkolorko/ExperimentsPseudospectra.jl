using ExperimentsPseudospectra
using Documenter

DocMeta.setdocmeta!(ExperimentsPseudospectra, :DocTestSetup, :(using ExperimentsPseudospectra); recursive=true)

makedocs(;
    modules=[ExperimentsPseudospectra],
    authors="Isaia Nisoli",
    sitename="ExperimentsPseudospectra.jl",
    format=Documenter.HTML(;
        canonical="https://orkolorko.github.io/ExperimentsPseudospectra.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/orkolorko/ExperimentsPseudospectra.jl",
    devbranch="main",
)
