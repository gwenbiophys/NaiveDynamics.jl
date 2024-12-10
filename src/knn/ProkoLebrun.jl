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
    build_bvh,
    traverse_bvh,
    batch_build_traverse,
    batched_batch_build

# how to build towards an API that makes it easy to extend a BVH procedure to different kinds 
# of shapes and considering different distances, as in the Noneuclidean paper?


struct SpheresBVHSpecs{T, K} <: SimulationSpecification
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


function SpheresBVHSpecs(; floattype, interaction_distance, leaves_count )
    if leaves_count < 2
        error("Please use more than one leaf")
    end
    branches_count = leaves_count - 1

    # this is arbitrary and really just my personal demonstration of struct instantiation with same name functions
    if floattype==Float32
        morton_type = Int32
        morton_length = Int32(10)

        return SpheresBVHSpecs{floattype, morton_type}(interaction_distance, leaves_count, branches_count, morton_length)
    elseif floattype==Float64
        morton_type = Int64
        morton_length = Int64(21)

        return SpheresBVHSpecs{floattype, morton_type}(interaction_distance, leaves_count, branches_count, morton_length)
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
struct QuantizedAABB{T} <: AxisAlignedBoundingBox
    index::T
    centroid::MVector{3, T}
    min::MVector{3, T}
    max::MVector{3, T}
end

#TODO does this need to be mutable?
# can index of the original atom be stored in a companion array? would this matter?
mutable struct GridKey{T, K} <: AABBGridKey
    index::T
    morton_code::T
    min::MVector{3, K}
    max::MVector{3, K}
    left::T
    skip::T
end

function generate_aabb(position::Vec3D{T}, spec::SpheresBVHSpecs{T}) where T
    aabb_array = [AABB{T}(i, position[i], position[i] .- spec.critical_distance, position[i] .+ spec.critical_distance) for i in eachindex(position)]


    return aabb_array

end
# this should be passed directly to the leaves.
function update_aabb!(position::Vec3D{T}, spec::SpheresBVHSpecs{T}, aabb_array::Vector{AABB{T}}) where T
    for i in eachindex(aabb_array)
        aabb_array[i].centroid .= position[i]
        aabb_array[i].min .= position[i] .- spec.critical_distance
        aabb_array[i].max .= position[i] .+ spec.critical_distance
    end
    
    return aabb_array
end


function update_gridkeys!(quantized_aabbarray, aabb_array, spec::SpheresBVHSpecs{T}, morton_type, clct) where T

    for dim in eachindex(clct.minDim)
        sort!(aabb_array, by = x -> x.min[dim])
        for i in eachindex(aabb_array)
            quantized_aabbarray[aabb_array[i].index].min[dim] = i 
        end

    end

    for dim in eachindex(clct.maxDim)
        sort!(aabb_array, by = x -> x.max[dim])
        for i in eachindex(aabb_array)
            quantized_aabbarray[aabb_array[i].index].max[dim] = i
        end
    end

    for dim in eachindex(clct.maxDim)
        sort!(aabb_array, by = x -> x.centroid[dim])
        for i in eachindex(aabb_array)
            quantized_aabbarray[aabb_array[i].index].centroid[dim] = i
        end
    end


end

"""
    update_mortoncodes!(L, quantized_aabbs, morton_length, morton_type)

Take an array of GridKeys, L, an array of 3D integer coordinates, quantized aabbs, and specification information,
to generate morton codes for each GridKey.
"""
# TODO how can we dump all of the correct digits from a source number into the product number
# without having to sequentially shift the bits back and forth
function update_mortoncodes!(L, quantized_aabbs, morton_length, morton_type) 
    # TODO should these be folded into the spec?

    inbit = morton_type(0)
    t3 = morton_type(3)
    t1 = morton_type(1)
    #L is an array of grid keys with an 'index' field which points to the 'nth index of a vector in the objectCollection struct' 
        # or a particular atom
    #n is the current nth bit of our morton code to change we wish to change, and it corresponds with every 3rd bit of our grid positions
    Threads.@threads for each in eachindex(quantized_aabbs)
        # set morton code to zero allows for data reuse
        L[each].morton_code = morton_type(0)

        for m in morton_type(morton_length):-1:t1 #iterate backwards
            for dim in t3:-1:t1

                # shift left to 'cancel' out the values of all bits before the m'th bit
                # then shift all the way to the right, filling in zeros
                inbit = (quantized_aabbs[each].centroid[dim] << (32 - m)) >>> 31

                # shift the morton code left by 1 bit to prevent overwrite
                # if inbit is 1, then morton_code becomes 1 at the right end
                # if inbit is 0, then do nothing, morton_code is already correct
                L[each].morton_code = (L[each].morton_code << 1) | inbit

            end
        end
    end
end




function reverse_bit(n::Int32)
    ret, power = 0, 31
    while n != 0
        ret += (n & 1) << power
        power -= 1
        n = n >> 1
    end

    return ret
end
function sort_mortoncodes!(L::Vector{GridKey{T, K}}) where {T, K}
    #make a specialized radix sort to replace base sort // space for GPU backends to put forward their own float
    sort!(L, by = x -> bitstring(x.morton_code), rev=false) # sorts lexicographically both the binary and the integer
    #sort!(L, by=x -> count(c -> c == '1', bitstring(x.morton_code)))
    #sort!(L, by = x -> x.morton_code), alg=RadixSort #wont run, RadixSort does not have iterate defined
end

function create_mortoncodes(position, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T}) where {T, K}
    #TODO is the kind of function scoping we want? -- ask again in the refactor
    aabb_array = generate_aabb(position, spec)
    morton_length = 0::Int64
    morton_string = " "::String
    #morton_type = K
    #idk how this should be better done to prevent runtime evaluation of a stupid if statement


    quantized_aabbs = [QuantizedAABB{K}(i, MVector{3, K}(0, 0, 0),  MVector{3, K}(0, 0, 0), MVector{3, K}(0, 0, 0)) for i in 1:spec.leaves_count]::Vector{QuantizedAABB{K}}

    update_gridkeys!(quantized_aabbs, aabb_array, spec, K, clct)

    L = [GridKey{K, T}(quantized_aabbs[i].index, 0, aabb_array[i].min, aabb_array[i].max, 0, 0) for i in 1:spec.leaves_count]

    update_mortoncodes!(L, quantized_aabbs, spec.morton_length, K)

    sort_mortoncodes!(L)

    return L

end


###### Phase 2: tree construction


function delta(i, L::Vector{GridKey{K, T}}, spec::SpheresBVHSpecs{T,K}) where {T, K}

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


#TODO what is the best way to do this? This is too naive.

#     expand_volume!(dst, operator, operand_a, operand_b) -> dst
# Conditionally expand the volume of x, y, z components of 'dst', where dst is expected to be either operand_a or b.

function expand_volume!(dst, operator, operand_a, operand_b)
    #copyto!(bounding_volume, expander)
    for dim in eachindex(bounding_volume)
        dst[dim] = operator(operand_a[dim], operand_b[dim]) * expander[dim]
    end


end

function stackless_interior!(store::Vector{Base.Threads.Atomic{Int64}}, i, nL, nI, keys, spec::SpheresBVHSpecs{T, K}) where {T, K}
    rangel = i 
    ranger = i
    dell = delta(rangel - 1, keys, spec)
    delr = delta(ranger, keys, spec)
    bounding_volume = [keys[i].min, keys[i].max]




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

            for dim in eachindex(bounding_volume[1])
                bounding_volume[2][dim] = (bounding_volume[2][dim] < keys[rightChild].max[dim]) * keys[rightChild].max[dim]
            end


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


            for dim in eachindex(bounding_volume[1])
                bounding_volume[1][dim] = (bounding_volume[1][dim] > keys[leftChild].min[dim]) * keys[leftChild].min[dim]
                #operator(operand_a[dim], operand_b[dim]) * expander[dim]
            end

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


function update_stackless_bvh!(keys, spec::SpheresBVHSpecs{T, K}) where {T, K}

    # TODO this should be generated in build_bvh and reset in update_bvh!
    store = [Base.Threads.Atomic{Int64}(0) for i in 1:spec.branches_count]

    # TODO What is the best perf method of handling inlining and bounds checking here?
    # I don't want to macro my code to hell
    Threads.@threads for i in 1:spec.leaves_count #in perfect parallel
        stackless_interior!(store, i, spec.leaves_count, spec.branches_count, keys, spec)
    end

    #Naive method of resetting the route
    for dim in eachindex(keys[branch_index(1, spec)].min)
        keys[branch_index(1, spec)].min[dim] = T(0.0)
        keys[branch_index(1, spec)].max[dim] = T(1.0)
    end

end


###### Phase 3: traversal

function neighbor_traverse(keys, neighbors, position, spec::SpheresBVHSpecs{T, K}) where {T, K}
    currentKey = branch_index(1, spec)
    neighbors = [] #TODO best method for improving this?
    while true

        if currentKey == 0
            break
        end
    end
end


###### Phase 4: put it all together


function build_bvh(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T}) where {T, K}

    keys = create_mortoncodes(position, spec, clct)::Vector{GridKey{K, T}} 
    I = [GridKey{K, T}(0, 0, MVector{3, K}(0.0, 0.0, 0.0), MVector{3, K}(0.0, 0.0, 0.0), 0, 0) for i in 1:spec.branches_count]
    append!(keys, I)

    #bvh_solver!(L, I, spec)
    update_stackless_bvh!(keys, spec)

    return keys
    #return neighborlist = println(traverse_bvh1(position, L, I, spec))
    #return traverse_bvh1(position, L, I, spec)

end

#TODO this is just demonstrative for howw we have to fix this file's structure moving forward
# for the integer coordinate vector, we could pass on a tuple of refs to it and the gridkey array from build_bvh
function rebuild_bvh!(keys, position, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T}) where {T, K}
    aabb_array = generate_aabb(position, spec)
    morton_length = 0::Int64
    morton_string = " "::String
    #morton_type = K
    #idk how this should be better done to prevent runtime evaluation of a stupid if statement


    quantized_aabbs = [QuantizedAABB{K}(i, MVector{3, K}(0, 0, 0),  MVector{3, K}(0, 0, 0), MVector{3, K}(0, 0, 0)) for i in 1:spec.leaves_count]::Vector{QuantizedAABB{K}}

    update_gridkeys!(quantized_aabbs, aabb_array, spec, K, clct)

    update_mortoncodes!(keys, quantized_aabbs, spec.morton_length, K)

    sort_mortoncodes!(keys)

    update_stackless_bvh!(keys, spec)
end

function build_traverse_bvh(position, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T}) where {T, K}
    tree = build_bvh(position, spec, clct)
    neighbors = [] #TODO what's the best way to handle this initialization?
    neighbor_traverse(tree, neighbors, position, spec )
end