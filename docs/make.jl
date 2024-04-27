using Documenter
using NaiveDynamics

push!(LOAD_PATH,"../src/")
makedocs(sitename="MyAwesomePackage.jl Documentation",
         pages = [
            "Index" => "index.md",
            "An other page" => "anotherPage.md",
         ],
         format = Documenter.HTML(prettyurls = false)
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/gwenbiophys/NaiveDynamics.jl.git",
    devbranch = "main"
)
