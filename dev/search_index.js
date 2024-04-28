var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Order   = [:module, :type, :constant, :function, :macro]","category":"page"},{"location":"api/","page":"API","title":"API","text":"Modules = [NaiveDynamics]\nOrder   = [:module, :type, :constant, :function, :macro]","category":"page"},{"location":"api/#NaiveDynamics.Collector","page":"API","title":"NaiveDynamics.Collector","text":"Collector\n\nCollector super-type for simulation initialization.\n\n\n\n\n\n","category":"type"},{"location":"api/#NaiveDynamics.GenericRandomCollector","page":"API","title":"NaiveDynamics.GenericRandomCollector","text":"GenericRandomCollector(objectnumber, minsize, maxsize, minspeed, maxspeed)\n\nAn Collector-subtype meant to acquire additional information from the user about how to make their system. In the collection function, positions and velocities will be randomly seeded from this Collector's boundary values.\n\n\n\n\n\n","category":"type"},{"location":"api/#NaiveDynamics.collect_objects-Tuple{GenericRandomCollector}","page":"API","title":"NaiveDynamics.collect_objects","text":"collect_objects(Collector::GenericRandomCollector)\n\nReturn a GenericObjectCollection with positions and speeds randomly seeded, as specified by the Collector object\n\n\n\n\n\n","category":"method"},{"location":"devdiary/#Developer-Diary","page":"Developer Diary","title":"Developer Diary","text":"","category":"section"},{"location":"devdiary/#1.-Let's-document-27-April","page":"Developer Diary","title":"1. Let's document - 27 April","text":"","category":"section"},{"location":"devdiary/","page":"Developer Diary","title":"Developer Diary","text":"Today has been an interesting day. I have been working on my baseline code but really hardening my documentation so that I have a well collected space to start hucking my thoughts. ","category":"page"},{"location":"devdiary/","page":"Developer Diary","title":"Developer Diary","text":"For instance, I added CSV and NamedArrays .jl to prepare for testing on whether the ","category":"page"},{"location":"devdiary/","page":"Developer Diary","title":"Developer Diary","text":"mutable struct GenericObjectCollection <: ObjectCollection \n    name::AbstractArray{String, 1}\n\n    position::AbstractArray{AbstractFloat, 3}\n    velocity::AbstractArray{AbstractFloat, 3}\n\n    uniqueID::AbstractArray{UUID,1}\n\nend","category":"page"},{"location":"devdiary/","page":"Developer Diary","title":"Developer Diary","text":"makes any sense, or if I should place all this information into a single array. And then I could test how this version of GenericObjectCollection scales, comparing it against a NamedArray convention. The primary point is to minimize processing time on putting these arrays together, while keeping my function accesses to data meaningful. I do not want the following:","category":"page"},{"location":"devdiary/","page":"Developer Diary","title":"Developer Diary","text":"Collection = [AbstractArray{String, 1}, AbstractArray{AbstractFloat, 3}, AbstractArray{AbstractFloat, 3}]\nfunction do_something(Collection)\n    return Collection[1] - Collection[2]\nend","category":"page"},{"location":"devdiary/","page":"Developer Diary","title":"Developer Diary","text":"in which function writing depends on making sure I have the right index of my Collection, so I don't do something stupid, like add a force term to a velocity term. I admit, something like position = Collection[3] is reallly easy, but we are here for overengineered solutions. ;)","category":"page"},{"location":"devdiary/","page":"Developer Diary","title":"Developer Diary","text":"But then my julia package would not precompile! Because I added new packages without updating the project.toml, so I had to relearn how to add new dependencies to a package.","category":"page"},{"location":"devdiary/","page":"Developer Diary","title":"Developer Diary","text":"It has been an hour and I still haven't figured out how to automatically get docstrings for functions passed into an index/api page, like they have it over at Molly.jl. Oh well, it is likely far better to get the gh-pages version of the documentation working.","category":"page"},{"location":"","page":"Index","title":"Index","text":"../../README.md","category":"page"}]
}