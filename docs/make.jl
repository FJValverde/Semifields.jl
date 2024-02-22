using Documenter
using Semifields


push!(LOAD_PATH,"../src/")
makedocs(
    sitename = "Semifields.jl Documentation",
    pages = [
        "Index" => "index.md",
        "Another page" => "anotherPage.md",
    ],
    format = Documenter.HTML(prettyurls = false),
    modules = [Semifields]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/FJValverde/Semifields.jl.git",
    devbranch = "main"
)
