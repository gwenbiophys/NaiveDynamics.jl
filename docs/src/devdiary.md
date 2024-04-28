# Developer Diary

## 1. Let's document - 27 April
Today has been an interesting day. I have been working on my baseline code but really hardening my documentation so that I have a well collected space to start hucking my thoughts. 

For instance, I added CSV and NamedArrays .jl to prepare for testing on whether the 
```julia
mutable struct GenericObjectCollection <: ObjectCollection 
    name::AbstractArray{String, 1}

    position::AbstractArray{AbstractFloat, 3}
    velocity::AbstractArray{AbstractFloat, 3}

    uniqueID::AbstractArray{UUID,1}

end
```
makes any sense, or if I should place all this information into a single array. And then I could test how this version of GenericObjectCollection scales, comparing it against a NamedArray convention. The primary point is to minimize processing time on putting these arrays together, while keeping my function accesses to data meaningful. I do not want the following:

```julia
Collection = [AbstractArray{String, 1}, AbstractArray{AbstractFloat, 3}, AbstractArray{AbstractFloat, 3}]
function do_something(Collection)
    return Collection[1] - Collection[2]
end
```
in which function writing depends on making sure I have the right index of my Collection, so I don't do something stupid, like add a force term to a velocity term. I admit, something like ``` position = Collection[3]``` is reallly easy, but we are here for overengineered solutions. ;)


But then my julia package would not precompile! Because I added new packages without updating the project.toml, so I had to relearn how to add new dependencies to a package.



It has been an hour and I still haven't figured out how to automatically get docstrings for functions passed into an index/api page, like they have it over at Molly.jl. Oh well, it is likely far better to get the gh-pages version of the documentation working.