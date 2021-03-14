using CosmoMMF
using Documenter

DocMeta.setdocmeta!(CosmoMMF, :DocTestSetup, :(using CosmoMMF); recursive=true)

makedocs(;
    modules=[CosmoMMF],
    authors="James Sunseri",
    repo="https://github.com/James11222/CosmoMMF.jl/blob/{commit}{path}#{line}",
    sitename="CosmoMMF.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://James11222.github.io/CosmoMMF.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/James11222/CosmoMMF.jl",
)
