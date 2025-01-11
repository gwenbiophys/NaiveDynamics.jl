# NaiveDynamics Tutorial

## Bounding Volume Hierarchy Neighbor List

In order to generate a neighbor list using BVH traversal, we must
* aquire positional data
* initialize BVH data
* generate a BVH
* traverse the BVH

Here in Julia the full process is as follows:

```julia
f = jldopen("data/positions/positions.jld2", "r")
myposition = deepcopy(read(f, "pos5000"))
close(f)

bvhspec = SpheresBVHSpecs(; bounding_distance=0.2
                            threshold_distance=0.2, 
                            leaves_count=length(myposition),
                            floattype=Float32, 
)

treeData = TreeData(myposition, spec)

my_neighbor_list =  neighbor_traverse(treeData.tree, myposition, spec)
```

In the first three lines, we used JLD2.jl to open a file of pregenerated positions stored in the data folder of NaiveDynamics. We selected `pos5000` from the file, indicating a vector of 5000 mutable vectors, which are each 3 Float32's long. 

Next, we instantiated a specification that includes the `bounding_distance`. This value determines the bounding volume around each 'atom' or point in the `myposition` array. The `threshold_distance` is the maximum distance at which two points are considered neighbors. In this example, we leave these values the same but a simulation which rebuilds the BVH less than 1 time per step would have a `bounding_distance > threshold_distance`.

With a specification and an array of positions, we construct a bvh, along with related reusable data. 

Finally, neighbor traverse will test the overlap of each element of the `myposition` vector in a time efficient way, using the BVH constructed before, and return a list of neighbors. This list is a tuple with 3 components: index to first element of `myposition`, index to second element, and the Euclidean distance between them.