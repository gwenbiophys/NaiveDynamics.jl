###### Bottom-up/parallel construction of linear bounding volume hierarchies with stackless traversal
# as described by Andrey Prokopenko and Damien Lebrun-Gradie
# https://arxiv.org/abs/2402.00665
# and as implemented by them in ArborX
# https://github.com/arborx/ArborX


# For tree generation, I tried to follow ArborX's implementation as closely as possible.


# This implementation is focused on initialize once and update update update





export
    SpheresBVHSpecs,
    AxisAlignedBoundingBox,
    AABB,
    QuantizedAABB,
    GridKey,
    update_mortoncodes!,
    create_mortoncodes,
    branch_index,
    update_stackless_bvh!,
    neighbor_traverse,
    build_bvh,
    rebuild_bvh!,
    build_traverse_bvh




# how to build towards an API that makes it easy to extend a BVH procedure to different kinds 
# of shapes and considering different distances, as in the Noneuclidean paper?


struct SpheresBVHSpecs{T, K} #<: SimulationSpecification #broken due to load order in the main NaiveDynamics module
    critical_distance::T
    leaves_count::Int64
    branches_count::Int64
    morton_length::K

end
"""
    function SpheresBVHSpecs(; floattype, interaction_distance, atoms_count, bins_count=atoms_count )

Instantiate a specification towards a BVH of sphere primitives. The ```interaction_distance``` is the maximum interaction distance for the set of 
pairwise interactions that a BVH+traversal algorithm finds the neighbors of. The ```bins_count``` is the number of chunks each axis will be divided by
in order convert particles from real space to grid space. By default, ```bins_count = atoms count```. Morton_length. Though Howard et al., 2019 chose 1023 bins to fit each 
grid axis within a UInt8, instead of here where the integer that fits a grid space axis has the same number of bits as the ```floattype```. 
"""


function SpheresBVHSpecs(; floattype, critical_distance, leaves_count )
    if leaves_count < 2
        error("Please use more than one leaf")
    end
    branches_count = leaves_count - 1

    # this is arbitrary and really just my personal demonstration of struct instantiation with same name functions
    if floattype==Float32
        morton_type = Int32
        morton_length = Int32(10)

        return SpheresBVHSpecs{floattype, morton_type}(critical_distance, leaves_count, branches_count, morton_length)
    elseif floattype==Float64
        morton_type = Int64
        morton_length = Int64(21)

        return SpheresBVHSpecs{floattype, morton_type}(critical_distance, leaves_count, branches_count, morton_length)
    end
    
end

##### Phase 1: morton code construction and updating

#TODO this variable type needs to become either T or K (or something else) to be consistent with MDInput
abstract type AxisAlignedBoundingBox end 
struct AABB{T} <: AxisAlignedBoundingBox where T
    index::Int64
    centroid::MVector{3, T}
    min::MVector{3, T}
    max::MVector{3, T}
end

abstract type AABBGridKey end
# must be mutable. would probably be better if this included a mutable vector

struct QuantizedAABB{T} <: AxisAlignedBoundingBox
    index::T
    centroid::MVector{3, T}
    min::MVector{3, T}
    max::MVector{3, T}
end

#TODO does this need to be mutable? Can we improve perf by restructuring to be immutable
# can index of the original atom be stored in a companion array? would this matter?
mutable struct GridKey{T, K} <: AABBGridKey
    index::K
    morton_code::K
    min::MVector{3, T}
    max::MVector{3, T}
    left::K
    skip::K
end
struct XYZVectors{type}
    x::Vector{type}
    y::Vector{type}
    z::Vector{type}
end


function update_quantized_positions(index, quantized, perm, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T}) where {T, K}
    for each in eachindex(quantized.x)
        quantized.x[index.x[perm.x[each]]] = K(each) #TODO i hope this is right, i have no idea
    end
    for each in eachindex(quantized.y)
        quantized.y[index.y[perm.y[each]]] = K(each)
    end
    for each in eachindex(quantized.z)
        quantized.z[index.z[perm.z[each]]] = K(each)
    end

end
"""
    update_mortoncodes!(L, quantized_aabbs, morton_length, morton_type)

Take an array of GridKeys, L, an array of 3D integer coordinates, quantized aabbs, and specification information,
to generate morton codes for each GridKey.
"""

function update_mortoncodes!(L, quantized, perm, morton_length, morton_type) 
    inbit = morton_type(0)
    t3 = morton_type(3)
    t1 = morton_type(1)
    #L is an array of grid keys with an 'index' field which points to the 'nth index of a vector in the objectCollection struct' 
        # or a particular atom
    #n is the current nth bit of our morton code to change we wish to change, and it corresponds with every 3rd bit of our grid positions
    for each in eachindex(quantized.x)
        # set morton code to zero allows for data reuse
        L[each].morton_code = morton_type(0)

        for m in morton_type(morton_length):-1:t1 #iterate backwards
            #TODO this is deccelerated due to having to access 3 different arrays and having to address every bit of th morton codes individualy rather than as ensembles
            # in this case, quantized being an mvec or even svec of xyz dimensions would be helpful
            inbit = (quantized.x[each] << (32 - m)) >>> 31
            L[each].morton_code = (L[each].morton_code << 1) | inbit


            inbit = (quantized.y[each] << (32 - m)) >>> 31
            L[each].morton_code = (L[each].morton_code << 1) | inbit

            
            inbit = (quantized.z[each] << (32 - m)) >>> 31
            L[each].morton_code = (L[each].morton_code << 1) | inbit


        end
    end
end

function sort_mortoncodes!(L::Vector{GridKey{T, K}}, spec::SpheresBVHSpecs{T, K}) where {T, K}
    #make a specialized radix sort to replace base sort // space for GPU backends to put forward their own sort
    sort!(L, by = x -> x.morton_code, rev=false) # sorts lexicographically both the binary and the integer
    #sort!(L, by=x -> count(c -> c == '1', bitstring(x.morton_code)))
    #sort!(L, by = x -> x.morton_code), alg=RadixSort #wont run, RadixSort does not have iterate defined
end

function create_mortoncodes(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T}) where {T, K}
    position_xyz = XYZVectors{T}([position[i][1] for i in 1:spec.leaves_count],
                                 [position[i][2] for i in 1:spec.leaves_count],
                                 [position[i][3] for i in 1:spec.leaves_count]
    )
    index_xyz = XYZVectors{K}([i for i in 1:spec.leaves_count],
                              [i for i in 1:spec.leaves_count],
                              [i for i in 1:spec.leaves_count]
    )
    perm_xyz = XYZVectors{K}(
        sortperm(position_xyz.x),
        sortperm(position_xyz.y),
        sortperm(position_xyz.z)
    )
    #partialsort
    quantized_xyz = XYZVectors{K}(zeros(K, spec.leaves_count),
                              zeros(K, spec.leaves_count),
                              zeros(K, spec.leaves_count)
    )
    update_quantized_positions(index_xyz, quantized_xyz, perm_xyz, spec, clct)

    #TODO this has to be fixed with the permuters and syntax
    L = [GridKey{T, K}(i, 0, 
            #MVector{3, T}(position_xyz.x[perm_xyz.x[i]] - spec.critical_distance, position_xyz.y[perm_xyz.y[i]] - spec.critical_distance, position_xyz.z[perm_xyz.z[i]] - spec.critical_distance),
            #MVector{3, T}(position_xyz.x[perm_xyz.x[i]] + spec.critical_distance, position_xyz.y[perm_xyz.y[i]] + spec.critical_distance, position_xyz.z[perm_xyz.z[i]] + spec.critical_distance),
            position[i] .- spec.critical_distance, position[i] .+ spec.critical_distance,
            0, 0) for i in 1:spec.leaves_count
    ]

    update_mortoncodes!(L, quantized_xyz, perm_xyz, spec.morton_length, K)
    sort_mortoncodes!(L, spec)#@time "mortons" sort_mortoncodes!(L)
    store = [Base.Threads.Atomic{Int64}(0) for i in 1:spec.branches_count]
    
    return tuple(Ref(L), Ref(position_xyz), Ref(index_xyz), Ref(perm_xyz), Ref(quantized_xyz), Ref(store))
end


###### Phase 2: tree construction


function delta(i, L::Vector{GridKey{T, K}}, spec::SpheresBVHSpecs{T,K}) where {T, K}

    # maintain numinternal nodes as length(L)-1 for now bc i dont think itll make a difference
    if  i >= spec.leaves_count || i < 1 #|| j > length(L) || j < 1 # i dont know if the same operation should be done to i, but idk how else to fix
        return typemax(K)
    end

    # XOR the two numbers to get a number with bits set where a and b differ
    # transcribed from ArborX directly, treeconstruction
    x = xor(L[i].morton_code, L[i+1].morton_code)


    
    return x + (x == K(0)) * (typemin(K) + (K(i) ⊻ K(i+1))) - K(1)


end


function branch_index(a, spec::SpheresBVHSpecs{T, K}) where {T, K}
    return c = a + spec.leaves_count # this could do without the c = but i am not ready with testing yet
end


function stackless_interior!(store::Vector{Base.Threads.Atomic{Int64}}, i, nL, nI, keys, spec::SpheresBVHSpecs{T, K}) where {T, K}
    rangel = i 
    ranger = i
    dell = delta(rangel - 1, keys, spec)
    delr = delta(ranger, keys, spec)
    bounding_volume = [keys[i].min, keys[i].max] #TODO this should be doable without main memory allocations, but because MVec, must instantiate per leaf per stackless_interior! call






    # stackless lBVH requires that the 'right most' elements of the tree point to an artificial node
    if i == nL # 1-based indexing requires this to be nL, but 0-based would be number of internal nodes
        keys[i].skip = 0
    else
        # Does Leaf[i+1] have a more similar morton code to Leaf[i] or Leaf[i+2]?
        # If Leaf[i+1] and Leaf[i] are more similar,
        # Then Leaf[i+1] and Leaf[i] are the right and left children of the same internal node.
        ir = i + 1
        #TODO is the <= necessary, or is < good enough?
        if delr < delta(ir, keys, spec) # are i and i+1 siblings, or are i+1 and i+2 siblings?
            keys[i].skip = ir
        else
            
            keys[i].skip = branch_index(ir, spec)
            
        end
    end

    while true
        isLeftChild = delr < dell
        if isLeftChild
            leftChild = i
            
            split = ranger # split position between the range of keys covered by any given INode
            ranger = Threads.atomic_cas!(store[split], 0, rangel)
            copy!


            if ranger == 0 
                break # this is the first thread to have made it here, so it is culled
            end

            delr = delta(ranger, keys, spec)

            #here is wehre boundary computation is performed.
            # memory has to sync here for data safety, whatever that means

            rightChild = split + 1
            rightChildIsLeaf = (rightChild == ranger)
            #Threads.atomic_fence() # uncertain what this does and if it is necessary
            # oh this is going to be so damn gross
            if rightChildIsLeaf
                #expand_volume!(bounding_volume[2], keys[rightChild].max, >)
            else
                #expand_volume!(bounding_volume[2], keys[branch_index(rightChild, spec)].max, >)
                rightChild = branch_index(rightChild, spec)

            end


            for yep in eachindex(bounding_volume[2])
                if keys[rightChild].max[yep] < bounding_volume[2][yep]
                    bounding_volume[2][yep] = keys[rightChild].max[yep] 
                end
            end
            #bounding_volume[2] = (keys[rightChild].max .< bounding_volume[2]) .* keys[rightChild].max + (bounding_volume[2] .< keys[rightChild].max) .* bounding_volume[2]
            # for dim in eachindex(bounding_volume[1])
            #     bounding_volume[2][dim] = (bounding_volume[2][dim] < keys[rightChild].max[dim]) * keys[rightChild].max[dim]
            # end


        else

            split = rangel - 1
            rangel = Threads.atomic_cas!(store[split], 0, ranger) 


            if rangel == 0 
                break
            end

            dell = delta(rangel-1, keys, spec)

            leftChild = split

            leftChildIsLeaf = (leftChild == rangel)
            
            #unclear if this is necessary
            #Threads.atomic_fence()
            if leftChildIsLeaf
                #expand_volume!(bounding_volume[1], keys[leftChild].min, <)
            else
                leftChild = branch_index(leftChild, spec)
                #expand_volume!(bounding_volume[1], keys[leftChild].min, <)
            end



            #TODO can't this be done in a cooler way :(
            for yep in eachindex(bounding_volume[1])
                if keys[leftChild].min[yep] < bounding_volume[1][yep]
                    bounding_volume[1][yep] = keys[leftChild].min[yep] 
                end
            end
            #bounding_volume[1] = eval_triple .* keys[leftChild].min



            #bounding_volume[1] = (keys[leftChild].min .< bounding_volume[1]) .* keys[leftChild].min + (bounding_volume[1] .< keys[leftChild].min) .* bounding_volume[1]
            
            # for dim in eachindex(bounding_volume[1])
            #     bounding_volume[1][dim] = (bounding_volume[1][dim] > keys[leftChild].min[dim]) * keys[leftChild].min[dim]
            #     #operator(operand_a[dim], operand_b[dim]) * expander[dim]
            # end

        end


        q = delr < dell ? (ranger) : rangel # commented out because the expanded version is easier to debug
         
        parentNode = branch_index(q, spec)
        keys[parentNode].left = leftChild

        

        if ranger == nL
            keys[parentNode].skip = 0
        else
            r = ranger + 1
            if delr < delta(r, keys, spec)
                keys[parentNode].skip = r
            else
                keys[parentNode].skip = branch_index(r, spec) 
            end
        end

        # this is incomplete, as a parent node is not guaranteed to encapsulate the volumes of all children
        copyto!(keys[parentNode].min, bounding_volume[1])
        copyto!(keys[parentNode].max, bounding_volume[2])

        
        
        i = branch_index(q, spec)


        if i == branch_index(1, spec)
            return
        end
    end
end


function update_stackless_bvh!(keys, store, spec::SpheresBVHSpecs{T, K}) where {T, K}



    # TODO What is the best perf method of handling inlining and bounds checking here?
    # I don't want to macro my code to hell
    Threads.@threads for i in 1:spec.leaves_count #in perfect parallel
        stackless_interior!(store, i, spec.leaves_count, spec.branches_count, keys, spec)
    end

    #Naive method of resetting the root
    for dim in eachindex(keys[branch_index(1, spec)].min)
        keys[branch_index(1, spec)].min[dim] = T(0.0)
        keys[branch_index(1, spec)].max[dim] = T(1.0)
    end

    #reset the values of the store for next time
    fill!(store, Base.Threads.Atomic{Int64}(0))

end


###### Phase 3: traversal
function proximity_test!(neighbors::Vector{Tuple{K, K, T}},  query_index::K, currentKey::K, positions::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}
     
    if query_index == currentKey
        return
    else
        d2 = sqrt( sum((positions[currentKey] .- positions[query_index]) .^ 2))

        # only push new pairs that are close together and nonredundant
        if d2 <= spec.critical_distance && query_index < currentKey
            push!(neighbors, tuple(query_index, currentKey, d2))

        end
    end

end
function overlap_test(keys, currentKey, query_index, positions, spec::SpheresBVHSpecs{T, K}) where {T, K}
    return sum(keys[currentKey].min .< positions[query_index] .< keys[currentKey].max)
end

#TODO this can be made parallel by giving a neighbor chunk to each thread and mending these neighbor lists together at the end
function neighbor_traverse(keys::Vector{GridKey{T,K}}, positions::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}
    neighbors = Vector{Tuple{K, K, T}}(undef, 0)

    # clamp traversal to 1 index before the last leaf because if every other leaf has been considered, 
    # then the last leaf does not need to be reconsidered for a neighbor pair
    for query_index in range(start=K(1), stop=K(spec.branches_count))

        currentKey = branch_index(1, spec)

        while currentKey != 0 # currentKey is the sentinel, end traversal of the given query
            # does query at all overlap with the volume of currentKey
            overlap = overlap_test(keys, currentKey, query_index, positions, spec)




            if overlap > 0 
                if keys[currentKey].left == 0 # currentKey is a leaf node
                    proximity_test!(neighbors, query_index, currentKey, positions, spec)
                    currentKey = currentKey = keys[currentKey].skip
                else #currentKey is a branch node, traverse to the left
                    currentKey = keys[currentKey].left
                end
            else #query is not contained, can cut off traversal on the 'lefts' sequencef
                currentKey = keys[currentKey].skip
            end

        end
    end
    return neighbors
end


###### Phase 4: put it all together

function build_bvh(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T}) where {T, K}

    # store a tuple of references to each of the persistent data structures created during bvh construction
    # this feels more sane than tupling up a bunch arrays directly
    bvhData = create_mortoncodes(position, spec, clct)

    #keys = initializationData[1][]


    #TODO is this helpful? sometimes the initialization of I takes a very long time and 'z' should help
    I = [GridKey{T, K}(0, 0, MVector{3, T}(0.0, 0.0, 0.0), MVector{3, K}(0.0, 0.0, 0.0), 0, 0) for i in 1:spec.branches_count]

    append!(bvhData[1][], I)

    update_stackless_bvh!(bvhData[1][], bvhData[6][], spec)

    return bvhData
end

function rebuild_bvh!(treeData, position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T}) where {T, K}
    #update positions, convert from mutable vector to struct of vectors
    #::Tuple{Vector{GridKey{T, K}},XYZVectors{T}, XYZVectors{K}, XYZVectors{K}, XYZVectors{K} }
    #TODO can these loops be replaced?
    #fill!(treeData[2][].x, [position[each][1] for each in eachindex(position)])
    
    for each in eachindex(position)
        treeData[2][].x[each] = position[each][1]
        #copyto!(treeData[2][].x[each], position[each][1])
        #setindex!(treeData[2][].x, position[1])
    end
    for each in eachindex(position)
        treeData[2][].y[each] = position[each][2]
    end
    for each in eachindex(position)
        treeData[2][].z[each] = position[each][3]
    end
    # update permuters
    sortperm!(treeData[4][].x, treeData[2][].x)
    sortperm!(treeData[4][].y, treeData[2][].y)
    sortperm!(treeData[4][].z, treeData[2][].z)
    
    update_quantized_positions(treeData[3][], treeData[5][], treeData[4][], spec, clct)

    #update leaf boundaries based on new positions
    #this has worse allocation performance than the expanded version without syntactic sugar???
    for each in 1:spec.leaves_count
        treeData[1][][each].min .= position[each] .- spec.critical_distance

    end
    for each in 1:spec.leaves_count
        treeData[1][][each].max .= position[each] .+ spec.critical_distance #.= or = ?

    end
    # realign leaf boundaries with indices to the atom positions that they represent
    for each in 1:spec.leaves_count
        treeData[1][][each].index = each
    end

    update_mortoncodes!(treeData[1][], treeData[5][], treeData[4][], spec.morton_length, K)
    #sort_mortoncodes!(treeData[1][][1:spec.leaves_count], spec)#@time "mortons" sort_mortoncodes!(L)
    partialsort!(treeData[1][], 1:spec.leaves_count, by=x -> x.morton_code)
    
    update_stackless_bvh!(treeData[1][], treeData[6][], spec)
    
    #return tuple(Ref(L), Ref(position_xyz), Ref(index_xyz), Ref(perm_xyz), Ref(quantized_xyz))
    # neighbors = []

    #neighbor_traverse(treeData[1][],  position, spec)

    #return treeData
    
end

function build_traverse_bvh(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T}) where {T, K}
    treeData = build_bvh(position, spec, clct)
    keys = treeData[1][] 



    neighbors = neighbor_traverse(keys, position, spec)
        #neighbor_traverse(tree, neighbors, position[i], spec )

    return neighbors
end