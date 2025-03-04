# NaiveDynamics Tutorial

## Bounding Volume Hierarchy Neighbor List

In order to generate a neighbor list using BVH traversal, we must
* aquire positional data
* initialize BVH data
* generate a BVH
* traverse the BVH

Here in Julia the full process is as follows:

```julia
using NaiveDynamics
using JLD2

f = jldopen("data/positions/positions.jld2", "r")
myposition = deepcopy(read(f, "pos5000"))
close(f)

spec = SpheresBVHSpecs(; neighbor_distance=0.2, 
                            atom_count=length(myposition),
                            floattype=Float32,
                            atomsperleaf=4 
)


treeData = TreeData(myposition, spec)

my_neighbor_list =  neighbor_traverse(treeData.tree, treeData.position, spec)
```

In the first three lines, we used JLD2.jl to open a file of pregenerated positions stored in the data folder of NaiveDynamics. We selected `pos5000` from the file, indicating a vector of 5000 mutable vectors, which are each 3 Float32's long. 

The `neighbor_distance` is the maximum distance at which two points are considered neighbors. The `atomsperleaf` determines the number of morton-encoded and sorted atoms that will be assigned to each leaf. The `atomsperleaf' can be any whole number that divides evenly into the number of points, or atoms. Increasing `atomsperleaf` will reduce the height of the bvh, while increasing the amount of distance checks between points. The height of the bvh can be balanced with the `atom_count` specification to achieve optimal performance.

With a specification and an array of positions, we construct a bvh, along with related reusable data. 

Finally, neighbor traverse will test the overlap of each element of the `myposition` vector --which has been transformed into the position field of the ```TreeData``` structure-- in a time efficient way, using the BVH constructed before, and return a list of neighbors. This list is a tuple with 3 components: index to first element of `myposition`, index to second element, and the Euclidean distance between them.