###### Bottom-up/parallel construction of linear bounding volume hierarchies with stackless traversal
# as described by Andrey Prokopenko and Damien Lebrun-Gradie
# https://arxiv.org/abs/2402.00665
# and as implemented by them in ArborX
# https://github.com/arborx/ArborX


# For tree generation, I tried to follow ArborX's implementation as closely as possible.
# Beyond that, I've made it up as I went along.




export
    SpheresBVHSpecs,
    GridKey,
    update_mortoncodes!,
    TreeData,
    branch_index,
    update_stackless_bvh!,
    neighbor_traverse,
    expt_neighbor_traverse,
    build_bvh,
    rebuild_bvh!,
    build_traverse_bvh,
    IndexSafePosition


struct SpheresBVHSpecs{T, K} #<: SimulationSpecification #broken due to load order in the main NaiveDynamics module
    critical_distance::T
    leaves_count::Int64
    branches_count::Int64
    morton_length::K

end
"""
    SpheresBVHSpecs(; floattype, interaction_distance, atoms_count, bins_count=atoms_count )

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

#TODO does this need to be mutable? Can we improve perf by restructuring to be immutable
# can index of the original atom be stored in a companion array? would this matter?
mutable struct GridKey{T, K}
    index::K
    morton_code::K
    min::SVector{3, T}
    max::SVector{3, T}
    left::K
    skip::K
end
struct XYZVectors{type}
    x::Vector{type}
    y::Vector{type}
    z::Vector{type}
end
struct IndexSafePosition{T, K}
    index::K
    vec::SVector{3, T}
end

"""
    struct TreeData{T, K}
        tree::Vector{GridKey{T, K}}
        position::Vector{IndexSafePosition{T, K}}
        quantizedposition::XYZVectors{K}
        store::Vector{Base.Threads.Atomic{Int64}}
    end
All data necessary to create a bounding volume hierarchy for stackless traversal.
"""
struct TreeData{T, K}
    tree::Vector{GridKey{T, K}}
    position::Vector{IndexSafePosition{T, K}}
    quantizedposition::XYZVectors{K}
    store::Vector{Base.Threads.Atomic{Int64}}
end

"""
    update_mortoncodes!(L, quantized_aabbs, morton_length, morton_type)

Take an array of GridKeys, L, an array of 3D integer coordinates, quantized aabbs, and specification information,
to generate morton codes for each GridKey.
"""
function update_mortoncodes!(L, quantized, morton_length, morton_type) 
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




function update_quantized_positions_direct(pos, quantized, spec::SpheresBVHSpecs{T, K}) where {T, K}
    sort!(pos, by=x->x.vec[1])
    for each in eachindex(quantized.z)
        quantized.x[pos[each].index] = K(each)
    end
    sort!(pos, by=x->x.vec[2])
    for each in eachindex(quantized.y)
        quantized.y[pos[each].index] = K(each)
    end
    sort!(pos, by=x->x.vec[3])
    for each in eachindex(quantized.z)
        quantized.z[pos[each].index] = K(each)
    end
end


function TreeData(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}

    pos = [IndexSafePosition{T, K}(i, SVector{3, T}(position[i])) for i in eachindex(position)]

    quantized_xyz = XYZVectors{K}(zeros(K, spec.leaves_count),
                              zeros(K, spec.leaves_count),
                              zeros(K, spec.leaves_count)
    )
    update_quantized_positions_direct(pos, quantized_xyz, spec)


    L = [GridKey{T, K}(i, 0, 
            SVector{3, T}(position[i] .- spec.critical_distance), SVector{3, T}(position[i] .+ spec.critical_distance),
            0, 0) for i in 1:spec.leaves_count
    ]

    update_mortoncodes!(L, quantized_xyz, spec.morton_length, K)
    
    sort_mortoncodes!(L, spec)#@time "mortons" sort_mortoncodes!(L)
    store = [Base.Threads.Atomic{Int64}(0) for i in 1:spec.branches_count]

    return TreeData{T, K}(L, pos, quantized_xyz, store)
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
    return a + spec.leaves_count
end


function stackless_interior!(store::Vector{Base.Threads.Atomic{Int64}}, i, nL, nI, keys, spec::SpheresBVHSpecs{T, K}) where {T, K}
    rangel = i 
    ranger = i
    dell = delta(rangel - 1, keys, spec)
    delr = delta(ranger, keys, spec)
    boundingmin = keys[i].min
    boundingmax = keys[i].max





    # stackless lBVH requires that the 'right most' elements of the tree point to an artificial node
    if i == nL # 1-based indexing requires this to be nL, but 0-based would be number of internal nodes
        keys[i].skip = 0
    else
        # Does Leaf[i+1] have a more similar morton code to Leaf[i] or Leaf[i+2]?
        # If Leaf[i+1] and Leaf[i] are more similar,
        # Then Leaf[i+1] and Leaf[i] are the right and left children of the same internal node.
        ir = i + 1

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



            if ranger == 0 
                break # this is the first thread to have made it here, so it is culled
            end

            delr = delta(ranger, keys, spec)


            rightChild = split + 1
            rightChildIsLeaf = (rightChild == ranger)
            #Threads.atomic_fence() # uncertain what this does and if it is necessary

            if !rightChildIsLeaf
                rightChild = branch_index(rightChild, spec)
            end

            boundingmax = (keys[rightChild].max .< boundingmax) .* keys[rightChild].max .+ (boundingmax .< keys[rightChild].max) .* boundingmax

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
            if !leftChildIsLeaf
                leftChild = branch_index(leftChild, spec)
            end


            boundingmin = (keys[leftChild].min .< boundingmin) .* keys[leftChild].min .+ (boundingmin .< keys[leftChild].min) .* boundingmin


        end


        q = delr < dell ? (ranger) : rangel 
         
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

        keys[parentNode].min = boundingmin
        keys[parentNode].max = boundingmax

        
        
        i = branch_index(q, spec)


        if i == branch_index(1, spec)
            return
        end
    end
end


function update_stackless_bvh!(keys, store, spec::SpheresBVHSpecs{T, K}) where {T, K}

    Threads.@threads for i in 1:spec.leaves_count #in perfect parallel
        stackless_interior!(store, i, spec.leaves_count, spec.branches_count, keys, spec)
    end

    keys[branch_index(1, spec)].min = SVector{3, T}(0.0, 0.0, 0.0)
    keys[branch_index(1, spec)].max = SVector{3, T}(1.0, 1.0, 1.0)

end


###### Phase 3: traversal
@inline function proximity_test!(neighbors::Vector{Tuple{K, K, T}},  query_index::K, currentKey::K, positions::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}
    
    # eliminate redundant pairings when qi == cK 
    #AND when [qi = a, cK = b] with [qi = b, cK = a] appear in the same pairs list array
    if !(query_index < currentKey) 
        return
    else
        d2 = sqrt( sum((positions[currentKey] .- positions[query_index]) .^ 2))

        # only push! new pairs that are close together 
        if d2 <= spec.critical_distance
            push!(neighbors, tuple(query_index, currentKey, d2))

        end
    end

end
@inline function overlap_test(keys, currentKey, query_index, positions, spec::SpheresBVHSpecs{T, K}) where {T, K}
    return sum(keys[currentKey].min .< positions[query_index] .< keys[currentKey].max)
end

@inline function newoverlap_test(keys, currentKey, query_index, positions, spec::SpheresBVHSpecs{T, K}) where {T, K}
    return sum(keys[currentKey].min .< positions[query_index] .< keys[currentKey].max) 
end


function neighbor_traverse(keys::Vector{GridKey{T,K}}, positions::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}

    threads = K(Threads.nthreads())

    #TODO isn't this structure extremely unfriendly to growing the individual elements at different times to different extents?
    # ----- not sure, having this as a vector of references seemed to diminish performance
    neighbor_vec = [Vector{Tuple{K, K, T}}(undef, 0) for i in 1:threads]


    # clamp traversal to 1 index before the last leaf because if every other leaf has been considered, 
    # then the last leaf does not need to be reconsidered for a neighbor pair
    Threads.@threads for chunk in 1:threads
        for query_index in K(chunk):threads:K(length(positions)-1)#range(start=K(1), stop=K(spec.branches_count))
            currentKey = branch_index(1, spec)

            while currentKey != 0 # currentKey is the sentinel, end traversal of the given query
                # does query at all overlap with the volume of currentKey
                overlap = overlap_test(keys, currentKey, query_index, positions, spec)





                if overlap > 0 
                    if keys[currentKey].left == 0 # currentKey is a leaf node
                        proximity_test!(neighbor_vec[chunk], query_index, currentKey, positions, spec)
                        currentKey = currentKey = keys[currentKey].skip
                    else #currentKey is a branch node, traverse to the left
                        currentKey = keys[currentKey].left
                    end
                else #query is not contained, can cut off traversal on the 'lefts' sequencef
                    currentKey = keys[currentKey].skip
                end

            end
        end
    end
    neighbors = neighbor_vec[1]
    for each in 2:1:threads
        append!(neighbors, neighbor_vec[each] )
    end

    return neighbors
end

function expt_neighbor_traverse(keys::Vector{GridKey{T,K}}, positions::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}

    threads = K(Threads.nthreads())

    #TODO isn't this structure extremely unfriendly to growing the individual elements at different times to different extents?
    # ----- not sure, having this as a vector of references seemed to diminish performance
    neighbor_vec = [Vector{Tuple{K, K, T}}(undef, 0) for i in 1:threads]


    # clamp traversal to 1 index before the last leaf because if every other leaf has been considered, 
    # then the last leaf does not need to be reconsidered for a neighbor pair
#     currentKey = 5
#     query_index = 2

#     overlap = @btime overlap_test($keys, $currentKey, $query_index, $positions, $spec)
#    # overlapb = @btime newoverlap_test($keys, $currentKey, $query_index, $positions, $spec)
#     println()
#     return

    Threads.@threads for chunk in 1:threads
        for query_index in K(chunk):threads:K(length(positions)-1)#range(start=K(1), stop=K(spec.branches_count))
            currentKey = branch_index(1, spec)

            while currentKey != 0 # currentKey is the sentinel, end traversal of the given query
                # does query at all overlap with the volume of currentKey
                overlap = newoverlap_test(keys, currentKey, query_index, positions, spec)




                if overlap > 0 
                    if keys[currentKey].left == 0 # currentKey is a leaf node
                        proximity_test!(neighbor_vec[chunk], query_index, currentKey, positions, spec)
                        currentKey = currentKey = keys[currentKey].skip
                    else #currentKey is a branch node, traverse to the left
                        currentKey = keys[currentKey].left
                    end
                else #query is not contained, can cut off traversal on the 'lefts' sequencef
                    currentKey = keys[currentKey].skip
                end

            end
        end
    end
    neighbors = neighbor_vec[1]
    for each in 2:1:threads
        append!(neighbors, neighbor_vec[each] )
    end

    return neighbors
end

###### Phase 4: put it all together

function build_bvh(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}

    treeData = TreeData(position, spec)

    I = [GridKey{T, K}(0, 0, SVector{3, T}(0.0, 0.0, 0.0), SVector{3, K}(0.0, 0.0, 0.0), 0, 0) for i in 1:spec.branches_count]
    
    append!(treeData.tree, I)

    update_stackless_bvh!(treeData.tree, treeData.store, spec)

    return treeData
end

function rebuild_bvh!(treeData, position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}

    for i in eachindex(position)
        treeData.position[i] = IndexSafePosition{T,K}(i, SVector{3, T}(position[i]))
    end

    update_quantized_positions_direct(treeData.position, treeData.quantizedposition, spec)


    #update leaf boundaries based on new positions
    #this has worse allocation performance than the expanded version without syntactic sugar???
    for each in 1:spec.leaves_count
        treeData.tree[each].min = SVector{3, T}(position[each] .- spec.critical_distance)

    end
    for each in 1:spec.leaves_count
        treeData.tree[each].max = SVector{3, T}(position[each] .+ spec.critical_distance) #.= or = ?

    end
    # realign leaf boundaries with indices to the atom positions that they represent
    for each in 1:spec.leaves_count
        treeData.tree[each].index = each
    end
    
    # have to reset to zero because ( I believe) zero values are not set
    # instead are unchanged from initialization
    # thus, we have to reset here
    for each in eachindex(treeData.tree) 
        treeData.tree[each].left = 0
        treeData.tree[each].skip = 0
    end

    update_mortoncodes!(treeData.tree, treeData.quantizedposition, spec.morton_length, K)

    leaves = treeData.tree[1:spec.leaves_count]
    sort!(leaves, by = x -> x.morton_code)
    treeData.tree[1:spec.leaves_count] = leaves #lmfao



    #reset boundaries
    for each in spec.leaves_count+1:1:spec.leaves_count+spec.branches_count
        treeData.tree[each].min = SVector{3, T}(0.0, 0.0, 0.0)
        treeData.tree[each].max = SVector{3, T}(0.0, 0.0, 0.0)
    end

    #reset store, have to dereference as they are atomic values
    for each in eachindex(treeData.store)
        treeData.store[each][] = 0
    end

    update_stackless_bvh!(treeData.tree, treeData.store, spec)
    
end

function build_traverse_bvh(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}
    treeData = build_bvh(position, spec)
    return neighbor_traverse(treeData.tree, position, spec)
end