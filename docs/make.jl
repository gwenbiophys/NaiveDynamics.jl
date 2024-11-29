
using Documenter
using NaiveDynamics

#push!(LOAD_PATH,"../src/")
makedocs(sitename="NaiveDynamics.jl",
         pages = [
            "Home" => "index.md",
            "Let's Talk Bounding Volume Hierarchies" => "bvh.md",
            "Developer Diary" => "devdiary.md",
            "Development Directions" => "devdirections.md",
            "Help Gwen" => "helpGwen.md"
         ],
         format = Documenter.HTML(prettyurls=true),

)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/gwenbiophys/NaiveDynamics.jl.git"
)

