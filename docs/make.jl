push!(LOAD_PATH,"C:/Users/Trist/.julia/dev/NaiveMD/NaiveDynamics.jl/src/")
using Documenter
using NaiveDynamics

#push!(LOAD_PATH,"../src/")
makedocs(sitename="NaiveDynamics.jl Documentation",
         pages = [
            "Index" => "index.md",
         ],
         format = Documenter.HTML(prettyurls = true)
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
"""
deploydocs(
    repo = "github.com/gwenbiophys/NaiveDynamics.jl.git",
    devbranch = "gh-pages"
)
"""
