###### Bottom-up/parallel construction of linear bounding volume hierarchies with stackless traversal
# as described by Andrey Prokopenko and Damien Lebrun-Gradie
# https://arxiv.org/abs/2402.00665
# and as implemented by them in ArborX
# https://github.com/arborx/ArborX


# For tree generation, I tried to follow ArborX's implementation as closely as possible.
# Beyond that, I've made it up as I went along.

# GridKey compression into bin coordinates was implemented as described in Howard et al. 2019
#https://doi.org/10.1016/j.commatsci.2019.04.004




export
    SpheresBVHSpecs,
    exptGridKey,
    GridKey,
    exptmortoncodes!,
    mortoncodes!,
    exptTreeData,
    TreeData,
    TreeData!,
    branch_index,
    bounding_volume_hierarchy!,
    neighbor_traverse,
    expt_neighbor_traverse,
    build_traverse_bvh,
    exptbuild_traverse_bvh,
    mortoncodes!,
    IPointPrimitive,
    proximity_test!,
    exptquantized_positions!,
    sort_mortoncodes!,
    exptcluster_primitives,
    exptbounding_volume_hierarchy!,
    binwidth, 
    positionToInt,
    intToPosition

"""
struct SpheresBVHSpecs{T, K} 
    neighbor_distance::T
    atom_count::Int64
    leaves_count::Int64
    branches_count::Int64
    atomsperleaf::Int64


Instantiate a specification towards a BVH of sphere primitives. 
The ```neighbor_distance``` is the distance at which two point primitives will be evaluated as neighbors when generating a neighbor list, using the 
```proximity_test!```` function. The ```leaves_count``` is the number of leaves, determined as the ratio of ```atom_count``` and ```atomsperleaf```.
"""
struct SpheresBVHSpecs{T, K} 
    neighbor_distance::T
    atom_count::Int64
    leaves_count::Int64
    branches_count::Int64
    atomsperleaf::Int64
end

function SpheresBVHSpecs(; neighbor_distance, atom_count, floattype, atomsperleaf )


    if rem(atom_count, atomsperleaf) != 0
        error("Please use an 'atomsperleaf' that evenly divides into 'atom_count' in BVH Specification")
    end
    leaves_count = atom_count / atomsperleaf
    if leaves_count < 2
        error("Please use more than one leaf in BVH Specification")
    end
    branches_count = leaves_count - 1

    # this is arbitrary and really just my personal demonstration of struct instantiation with same name functions
    if floattype==Float32
        morton_type = Int32

        return SpheresBVHSpecs{floattype, morton_type}(neighbor_distance, atom_count, leaves_count, branches_count, atomsperleaf)
    elseif floattype==Float64
        morton_type = Int64

        return SpheresBVHSpecs{floattype, morton_type}(neighbor_distance, atom_count, leaves_count, branches_count, atomsperleaf)
    end
    
end



##### Phase 1: morton code construction and updating

#TODO does this need to be mutable? Can we improve perf by restructuring to be immutable
# can index of the original atom be stored in a companion array? would this matter?
struct GridKey{T, K}
    #index::K
    #morton_code::K
    min::SVector{3, T}
    max::SVector{3, T}
    left::K
    skip::K
end
struct exptGridKey{T, K}
    min::Int32 #in a float64 scenario, Int64 would not capture anywhere near the ratio of significance that Int32 can capture from 3 float32s
    max::Int32
    left::K
    skip::K
end

#TODO deprecated in name of other methods
struct XYZVectors{type}
    x::Vector{type}
    y::Vector{type}
    z::Vector{type}
end
struct IndexSafePosition{T, K}
    index::K
    vec::SVector{3, T}
end

# eachh point primitive has an 'index' to its original index place in the position (force and velocity) array(s)
##TODO expt if this can be immutable
# mutable struct PointPrimitive{T, K}
#     index::K
#     morton_code::K #TODO update naming so this becomes one word, maybe even just morton?
#     position::SVector{3, T}
# end
struct IPointPrimitive{T, K}
    index::K
    morton_code::K #TODO update naming so this becomes one word, maybe even just morton?
    position::SVector{3, T}
end

#TODO deprecated in the name of tuple
"""
    struct TreeData{T, K}
        tree::Vector{GridKey{T, K}}
        position::Vector{IndexSafePosition{T, K}}
        quantizedposition::Vector{MVector{3, K}}
        store::Vector{Base.Threads.Atomic{Int64}}
    end
Generate a stackless bounding volume hierarchy stored in the `tree` field, while storing related data in tree construction for later reuse.
"""
# struct TreeData{T, K}
#     tree::Vector{MGridKey{T, K}}
#     position::Vector{PointPrimitive{T, K}}
#     quantizedposition::Vector{MVector{3, K}}
#     store::Vector{Base.Threads.Atomic{K}}
# end
# struct exptTreeData{T, K}
#     tree::Vector{GridKey{T, K}}
#     position::Vector{IPointPrimitive{T, K}}
#     quantizedposition::Vector{MVector{3, K}}
#     store::Vector{Base.Threads.Atomic{K}}

# end
# struct newexptTreeData{T, K}
#     tree::Ref{Vector{GridKey{T, K}}}
#     position::Ref{Vector{IPointPrimitive{T, K}}}
#     quantizedposition::Ref{Vector{MVector{3, K}}}
#     store::Ref{Vector{Base.Threads.Atomic{K}}}

# end

const binwidth = Float32(1/1023)

#TODO i dont like these names, they are generic, but also too long?
#but also shrink and grow dont feel right
function compress_position(pos::SVector{3, Float32}, UpOrDown)
    intermediate = round.(Int32, pos ./ binwidth, UpOrDown)
    result = Int32(0)
    for i in eachindex(intermediate)
        result = result << Int32(10) | intermediate[i]

    end

    return result
end

function expand_integer(num::Int32)
    return SVector{3, Float32}( (num << Int32(2+10*(i-1))) >>> Int32(22)  * binwidth for i in 1:3)
end


"""
    mortoncodes!(L, quantized_aabbs, spec::SpheresBVHSpecs{T, K}) where {T, K} )

Take an array of GridKeys, L, an array of 3D integer coordinates, quantized aabbs, and specification information,
to generate morton codes for each GridKey.

For a discussion into how this function works, see 11 January, 2025 in devdiary, approx. L1245.
"""

#TODO remove L notation, is confusing and incorrect
function exptmortoncodes!(L, quantized, spec::SpheresBVHSpecs{T, K}) where {T, K} 
    #TODO implement for morton_type::Int32, and Int64 as the magic numbers change
    #TODO is it possible to avoid the magic not values? so that this runs on selective promotion of zero to 1
    # i.e. a morton code bit will be left as zero because it was not modified, but becomes 1 in order to fulfill bit interleaving?
    magic_values = SVector{3, Int32}(153391689, 306783378, 613566756)
    magic_not_values = SVector{3, Int32}(-153391690, -306783379, -613566757 )

    for each in eachindex(quantized)
        # set morton code to zero allows for data reuse
        input= K(0)

        for dim in eachindex(quantized[1])
            yin = quantized[each][dim] & magic_values[dim]
            yang = quantized[each][dim] | magic_not_values[dim]
            yinyang = yin & yang
            input |= yinyang

        end

        L[each] = IPointPrimitive{T,K}(L[each].index, input, L[each].position)
    end

end
function mortoncodes!(L, quantized, spec::SpheresBVHSpecs{T, K}) where {T, K} 
    #TODO implement for morton_type::Int32, and Int64 as the magic numbers change
    #TODO is it possible to avoid the magic not values? so that this runs on selective promotion of zero to 1
    # i.e. a morton code bit will be left as zero because it was not modified, but becomes 1 in order to fulfill bit interleaving?
    magic_values = SVector{3, Int32}(153391689, 306783378, 613566756)
    magic_not_values = SVector{3, Int32}(-153391690, -306783379, -613566757 )

    for each in eachindex(quantized)
        # set morton code to zero allows for data reuse
        input= K(0)

        for dim in eachindex(quantized[1])
            yin = quantized[each][dim] & magic_values[dim]
            yang = quantized[each][dim] | magic_not_values[dim]
            yinyang = yin & yang
            input |= yinyang

        end

        L[each] = IPointPrimitive{T,K}(L[each].index, input, L[each].position)
    end
end
function sort_mortoncodes!(L, spec::SpheresBVHSpecs{T, K}) where {T, K}
    #make a specialized radix sort to replace base sort // space for GPU backends to put forward their own sort ;-)
    sort!(L, by = x -> x.morton_code, rev=false) # sorts lexicographically both the binary and the integer
    #sort!(L, by=x -> count(c -> c == '1', bitstring(x.morton_code)))
    #sort!(L, by = x -> x.morton_code), alg=RadixSort #wont run, RadixSort does not have iterate defined
end




function quantized_positions!(quantized::Vector{SVector{3, K}}, pos::Vector{IPointPrimitive{T,K}}, spec::SpheresBVHSpecs{T, K}) where {T, K}
    sort!(pos, by=x->x.position[1])
    for each in eachindex(pos)
        quantized[pos[each].index] = SVector{3, K}(K(each), quantized[pos[each].index][2], quantized[pos[each].index][3] )
    end

    sort!(pos, by=x->x.position[2])
    for each in eachindex(pos)
        quantized[pos[each].index] = SVector{3, K}(quantized[pos[each].index][1], K(each), quantized[pos[each].index][3] )
    end

    sort!(pos, by=x->x.position[3])
    for each in eachindex(pos)
        quantized[pos[each].index] = SVector{3, K}(quantized[pos[each].index][1],quantized[pos[each].index][2], K(each))
    end
end
function exptquantized_positions!(quantized::Vector{SVector{3, K}}, pos::Vector{IPointPrimitive{T,K}}, spec::SpheresBVHSpecs{T, K}) where {T, K}
    #TODO i think we can do better than have MVec, no?
    # for dim in eachindex(quantized[1])
    #     sort!(pos, by=x->x.position[dim])
    #     for each in eachindex(pos)
    #         quantized[pos[each].index][dim] = K(each)
    #     end
    # end

    sort!(pos, by=x->x.position[1])
    for each in eachindex(pos)
        quantized[pos[each].index] = SVector{3, K}(K(each), quantized[pos[each].index][2], quantized[pos[each].index][3] )
    end

    sort!(pos, by=x->x.position[2])
    for each in eachindex(pos)
        quantized[pos[each].index] = SVector{3, K}(quantized[pos[each].index][1], K(each), quantized[pos[each].index][3] )
    end

    sort!(pos, by=x->x.position[3])
    for each in eachindex(pos)
        quantized[pos[each].index] = SVector{3, K}(quantized[pos[each].index][1],quantized[pos[each].index][2], K(each))
    end
end


function cluster_primitives(L::Vector{IPointPrimitive{T, K}}, spec::SpheresBVHSpecs{T, K}) where {T, K}

    leaves = [GridKey{T,K}(SVector{3, T}(0.0, 0.0, 0.0), SVector{3, T}(0.0, 0.0, 0.0), 0, 0) for i in 1:spec.leaves_count]
    for each in eachindex(leaves)

        bird = (each-1) * spec.atomsperleaf + 1

        leaves[each] = GridKey{T,K}(#L[each].index, 
                                    #L[each].morton_code, 

                                    L[bird].position .- spec.neighbor_distance,
                                    L[bird].position .+ spec.neighbor_distance,
                                    0,
                                    0
        )

        for i in 1:1:spec.atomsperleaf

            a = (each-1)*spec.atomsperleaf  + i 


            leaves[each] = GridKey{T,K}(#L[each].index, 
                                    #L[each].morton_code, 

                                    min.(leaves[each].min, L[a].position .- spec.neighbor_distance), 
                                    max.(leaves[each].max, L[a].position .+ spec.neighbor_distance),
                                    0,
                                    0
            )
        end

    end

    return leaves
end

#TODO this function is entirely broken by compression method
function exptcluster_primitives(L::Vector{IPointPrimitive{T, K}}, spec::SpheresBVHSpecs{T, K}) where {T, K}

    leaves = [GridKey{T,K}(SVector{3, T}(0.0, 0.0, 0.0), SVector{3, T}(0.0, 0.0, 0.0), 0, 0) for i in 1:spec.leaves_count]
    for each in eachindex(leaves)

        bird = (each-1) * spec.atomsperleaf + 1

        leaves[each] = GridKey{T,K}(#L[each].index, 
                                    #L[each].morton_code, 
                                    L[bird].position,
                                    L[bird].position,
                                    # compress_position(L[bird].position , RoundDown),
                                    # compress_position(L[bird].position , RoundUp),
                                    0,
                                    0
        )

        for i in 1:1:spec.atomsperleaf

            a = (each-1)*spec.atomsperleaf  + i 

            #TODO it is probably worth breaking this out into separately evalulateable variables
            leaves[each] = GridKey{T,K}(#L[each].index, 
                                    #L[each].morton_code, 
                                    min.(leaves[each].min, L[a].position ), 
                                    max.(leaves[each].max, L[a].position),
                                    # compress_position(min.(expand_integer(leaves[each].min), L[a].position ), RoundDown), 
                                    # compress_position(max.(expand_integer(leaves[each].max), L[a].position ), RoundUp),
                                    0,
                                    0
            )
        end

    end

    return leaves
end

#TODO should TreeData be capitalized??? im waffling
function exptTreeData(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}


    pos = [IPointPrimitive{T,K}(i, 0, position[i]) for i in 1:spec.atom_count]

    quantized_xyz = [SVector{3, K}(0, 0, 0) for i in eachindex(position)] 

    exptquantized_positions!(quantized_xyz, pos, spec)

    exptmortoncodes!(pos, quantized_xyz, spec)

    sort_mortoncodes!(pos, spec)

    #1 alloc per item in this generator expression
    #store = [Threads.Atomic{K}(0) for i in 1:spec.branches_count]
    store = [K(0) for i in 1:spec.branches_count]
    # @time for each in eachindex(store)
    #     store[each] = Threads.Atomic{K}(store[each])
    # end


    L = exptcluster_primitives(pos, spec)



    I = [GridKey{T, K}(SVector{3, T}(0.0, 0.0, 0.0), SVector{3, K}(0.0, 0.0, 0.0), 0, 0) for i in 1:spec.branches_count]
    #I = [exptGridKey{T, K}(0, 0, 0, 0) for i in 1:spec.branches_count]
    
    append!(L, I)

    # i think this is the happiest median, particularly for CPU construction
    
    bounding_volume_hierarchy!(L, store, spec, pos)
    
    compressedL = [exptGridKey{T, K}(compress_position(L[i].min, RoundDown), compress_position(L[i].max, RoundUp), L[i].left, L[i].skip) for i in eachindex(L)]

    # for i in eachindex(L)
    #     println("min ", L[i].min)
    #     println(expand_integer(compressedL[i].min))

    #     println("max ", L[i].max)
    #     println(expand_integer(compressedL[i].max))
    #     println()
    # end

    
    return tuple(L, pos, quantized_xyz, store, compressedL)
end
function TreeData(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}

    pos = [IPointPrimitive{T,K}(i, 0, position[i]) for i in 1:spec.atom_count]

    quantized_xyz = [SVector{3, K}(0, 0, 0) for i in eachindex(position)] 

    quantized_positions!(quantized_xyz, pos, spec)

    mortoncodes!(pos, quantized_xyz, spec)

    sort_mortoncodes!(pos, spec)

    #1 alloc per item in this generator expression
    #store = [Threads.Atomic{K}(0) for i in 1:spec.branches_count]
    store = [K(0) for i in 1:spec.branches_count]
    # @time for each in eachindex(store)
    #     store[each] = Threads.Atomic{K}(store[each])
    # end


    L = cluster_primitives(pos, spec)




    I = [GridKey{T, K}(SVector{3, T}(0.0, 0.0, 0.0), SVector{3, K}(0.0, 0.0, 0.0), 0, 0) for i in 1:spec.branches_count]
    
    append!(L, I)


    bounding_volume_hierarchy!(L, store, spec, pos)


    return tuple(L, pos, quantized_xyz, store)
end



function TreeData!(treeData, position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}

    for i in eachindex(position)
        treeData.position[i] = PointPrimitive{T,K}(i, 0, SVector{3, T}(position[i]))
        #copyto!(treeData.position[i].vec, position[i])
    end

    quantized_positions!(treeData.quantizedposition, treeData.position, spec)


    #update leaf boundaries based on new positions
    #this has worse allocation performance than the expanded version without syntactic sugar???
    for each in 1:spec.leaves_count
        treeData.tree[each].min = SVector{3, T}(position[each] .- spec.neighbor_distance)

    end
    for each in 1:spec.leaves_count
        treeData.tree[each].max = SVector{3, T}(position[each] .+ spec.neighbor_distance) #.= or = ?

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

    mortoncodes!(treeData.tree, treeData.quantizedposition, spec)

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
        treeData.store[each][] = K(0)
    end

    #bounding_volume_hierarchy!(treeData.tree, treeData.store, spec)
    
end

###### Phase 2: tree construction

function delta(i, L, spec::SpheresBVHSpecs{T,K}) where {T, K}

    # maintain numinternal nodes as length(L)-1 for now bc i dont think itll make a difference
    if  i >= spec.leaves_count || i < 1 #|| j > length(L) || j < 1 # i dont know if the same operation should be done to i, but idk how else to fix
        return typemax(K)
    end

    # XOR the two numbers to get a number with bits set where a and b differ
    # transcribed from ArborX directly, treeconstruction
    #x = xor(L[i].morton_code, L[i+1].morton_code)
    atom_i = 1 + (i-1)*spec.atomsperleaf
    atom_i_and1 = 1 + (i)*spec.atomsperleaf
    x = xor(pos[atom_i].morton_code, pos[atom_i_and1].morton_code)


    #TODO is it supposed  to be i or ai here?    
    return x + (x == K(0)) * (typemin(K) + (K(i) ⊻ K(i+1))) - K(1)
end
function exptdelta(i, L, spec::SpheresBVHSpecs{T,K}, pos) where {T, K}

    # maintain numinternal nodes as length(L)-1 for now bc i dont think itll make a difference
    if  i >= spec.leaves_count || i < 1 #|| j > length(L) || j < 1 # i dont know if the same operation should be done to i, but idk how else to fix
        return typemax(K)
    end

    # XOR the two numbers to get a number with bits set where a and b differ
    # transcribed from ArborX directly, treeconstruction
    #x = xor(L[i].morton_code, L[i+1].morton_code)
    atom_i = 1 + (i-1)*spec.atomsperleaf
    atom_i_and1 = 1 + (i)*spec.atomsperleaf
    x = xor(pos[atom_i].morton_code, pos[atom_i_and1].morton_code)


    #TODO is it supposed  to be i or ai here?    
    return x + (x == K(0)) * (typemin(K) + (K(i) ⊻ K(i+1))) - K(1)


end


function branch_index(a, spec::SpheresBVHSpecs{T, K}) where {T, K}
    return K(a + spec.leaves_count)
end


Base.@propagate_inbounds function bvh_interior!(keys, store, i, nL, nI, spec::SpheresBVHSpecs{T, K}, pos) where {T, K}
    rangel = i 
    ranger = i
    dell = exptdelta(rangel - 1, keys, spec, pos)
    delr = exptdelta(ranger, keys, spec, pos)
    boundingmin = keys[i].min
    boundingmax = keys[i].max





    # stackless lBVH requires that the 'right most' elements of the tree point to an artificial node
    if i == nL # 1-based indexing requires this to be nL, but 0-based would be number of internal nodes
        #keys[i].skip = K(0)
        zerskip = K(0)
        #keys[i] = GridKey{T,K}(keys[i].index, keys[i].morton_code, keys[i].min, keys[i].max, keys[i].left, zerskip)
        keys[i] = GridKey{T,K}(keys[i].min, keys[i].max, keys[i].left, zerskip)
    else
        # Does Leaf[i+1] have a more similar morton code to Leaf[i] or Leaf[i+2]?
        # If Leaf[i+1] and Leaf[i] are more similar,
        # Then Leaf[i+1] and Leaf[i] are the right and left children of the same internal node.
        ir = i + K(1)

        if delr < exptdelta(ir, keys, spec, pos) # are i and i+1 siblings, or are i+1 and i+2 siblings?
            #keys[i].skip = ir
            zerskip = K(ir)
            #keys[i] = GridKey{T,K}(keys[i].index, keys[i].morton_code, keys[i].min, keys[i].max, keys[i].left, zerskip)
            keys[i] = GridKey{T,K}(keys[i].min, keys[i].max, keys[i].left, zerskip)
        else
            
            #keys[i].skip = branch_index(K(ir), spec)
            zerskip = branch_index(K(ir), spec)
            #keys[i] = GridKey{T,K}(keys[i].index, keys[i].morton_code, keys[i].min, keys[i].max, keys[i].left, zerskip)
            keys[i] = GridKey{T,K}(keys[i].min, keys[i].max, keys[i].left, zerskip)
            
        end
    end

    while true
        isLeftChild = delr < dell
        if isLeftChild
            leftChild = i

            split = ranger # split position between the range of keys covered by any given INode
            #ranger = Threads.atomic_cas!(store[split], K(0), rangel)
            rangereval = Atomix.@atomicreplace store[split] K(0) => rangel
            ranger = rangereval[1]

            if ranger == 0 
                break # this is the first thread to have made it here, so it is culled
            end

            delr = exptdelta(ranger, keys, spec, pos)


            rightChild = split + 1
            rightChildIsLeaf = (rightChild == ranger)
            #Threads.atomic_fence() # uncertain what this does and if it is necessary

            if !rightChildIsLeaf
                rightChild = branch_index(rightChild, spec)
            end


            boundingmax = max.(keys[rightChild].max, boundingmax)

            boundingmin = min.(keys[rightChild].min, boundingmin)
        else

            split = rangel - 1
            #rangel = Threads.atomic_cas!(store[split], K(0), ranger) 

            rangeleval = Atomix.@atomicreplace store[split] K(0) => ranger
            rangel = rangeleval[1]


            if rangel == 0 
                break
            end

            dell = exptdelta(rangel-1, keys, spec, pos)

            leftChild = split

            leftChildIsLeaf = (leftChild == rangel)
            
            #unclear if this is necessary
            #Threads.atomic_fence()
            if !leftChildIsLeaf
                leftChild = branch_index(leftChild, spec)
            end




            boundingmax = max.(keys[leftChild].max, boundingmax)

            boundingmin = min.(keys[leftChild].min, boundingmin)
        end


        q = delr < dell ? (ranger) : rangel 
         
        parentNode = branch_index(q, spec)
        #keys[parentNode].left = K(leftChild)
        pleftchild = leftChild
        #keys[parentNode] = GridKey{T,K}(keys[parentNode].index, keys[parentNode].morton_code, keys[parentNode].min, keys[parentNode].max, pleftchild, keys[parentNode].skip)
        keys[parentNode] = GridKey{T,K}(keys[parentNode].min, keys[parentNode].max, pleftchild, keys[parentNode].skip)
        

        if ranger == nL
            #keys[parentNode].skip = K(0)
            zerskip = K(0)
            #keys[parentNode] = GridKey{T,K}(keys[parentNode].index, keys[parentNode].morton_code, keys[parentNode].min, keys[parentNode].max, keys[parentNode].left, zerskip)
            keys[parentNode] = GridKey{T,K}(boundingmin, boundingmax, keys[parentNode].left, zerskip)
        else
            r = ranger + K(1)
            if delr < exptdelta(r, keys, spec, pos)
                #keys[parentNode].skip = r
                #keys[parentNode] = GridKey{T,K}(keys[parentNode].index, keys[parentNode].morton_code, keys[parentNode].min, keys[parentNode].max, keys[parentNode].left, r)
                keys[parentNode] = GridKey{T,K}(boundingmin, boundingmax, keys[parentNode].left, r)
            else
                #keys[parentNode].skip = branch_index(r, spec) 
                r = branch_index(r, spec) 
                #keys[parentNode] = GridKey{T,K}(keys[parentNode].index, keys[parentNode].morton_code, keys[parentNode].min, keys[parentNode].max, keys[parentNode].left, r)
                keys[parentNode] = GridKey{T,K}(boundingmin, boundingmax, keys[parentNode].left, r)
            end
        end

        # keys[parentNode].min = boundingmin
        # keys[parentNode].max = boundingmax

        #keys[parentNode] = GridKey{T,K}(keys[parentNode].index, keys[parentNode].morton_code, boundingmin, boundingmax, keys[parentNode].left, keys[parentNode].skip)
        #keys[parentNode] = GridKey{T,K}(boundingmin, boundingmax, keys[parentNode].left, keys[parentNode].skip)
        i = branch_index(q, spec)


        if i == branch_index(K(1), spec)
            return
        end
    end
end


function bounding_volume_hierarchy!(keys::Vector{GridKey{T,K}}, store, spec::SpheresBVHSpecs{T, K}, pos) where {T, K}

    @batch for i in 1:spec.leaves_count #in perfect parallel
        i = K(i)
        @inbounds bvh_interior!(keys, store, i, spec.leaves_count, spec.branches_count, spec, pos)
    end

    min = SVector{3, T}(0.0, 0.0, 0.0)
    max = SVector{3, T}(1.0, 1.0, 1.0)
    keys[branch_index(1, spec)] = GridKey{T,K}(min, max, keys[branch_index(1, spec)].left, keys[branch_index(1, spec)].skip)

end

Base.@propagate_inbounds function exptbvh_interior!(keys, store, i, nL, nI, spec::SpheresBVHSpecs{T, K}, pos) where {T, K}
    rangel = i 
    ranger = i
    dell = exptdelta(rangel - 1, keys, spec, pos)
    delr = exptdelta(ranger, keys, spec, pos)
    boundingmin = expand_integer(keys[i].min)
    boundingmax = expand_integer(keys[i].max)





    # stackless lBVH requires that the 'right most' elements of the tree point to an artificial node
    if i == nL # 1-based indexing requires this to be nL, but 0-based would be number of internal nodes
        #keys[i].skip = K(0)
        zerskip = K(0)
        #keys[i] = GridKey{T,K}(keys[i].index, keys[i].morton_code, keys[i].min, keys[i].max, keys[i].left, zerskip)
        keys[i] = GridKey{T,K}(keys[i].min, keys[i].max, keys[i].left, zerskip)
    else
        # Does Leaf[i+1] have a more similar morton code to Leaf[i] or Leaf[i+2]?
        # If Leaf[i+1] and Leaf[i] are more similar,
        # Then Leaf[i+1] and Leaf[i] are the right and left children of the same internal node.
        ir = i + K(1)

        if delr < exptdelta(ir, keys, spec, pos) # are i and i+1 siblings, or are i+1 and i+2 siblings?
            #keys[i].skip = ir
            zerskip = K(ir)
            #keys[i] = GridKey{T,K}(keys[i].index, keys[i].morton_code, keys[i].min, keys[i].max, keys[i].left, zerskip)
            keys[i] = GridKey{T,K}(keys[i].min, keys[i].max, keys[i].left, zerskip)
        else
            
            #keys[i].skip = branch_index(K(ir), spec)
            zerskip = branch_index(K(ir), spec)
            #keys[i] = GridKey{T,K}(keys[i].index, keys[i].morton_code, keys[i].min, keys[i].max, keys[i].left, zerskip)
            keys[i] = GridKey{T,K}(keys[i].min, keys[i].max, keys[i].left, zerskip)
            
        end
    end

    while true
        isLeftChild = delr < dell
        if isLeftChild
            leftChild = i

            split = ranger # split position between the range of keys covered by any given INode
            #ranger = Threads.atomic_cas!(store[split], K(0), rangel)
            rangereval = Atomix.@atomicreplace store[split] K(0) => rangel
            ranger = rangereval[1]

            if ranger == 0 
                break # this is the first thread to have made it here, so it is culled
            end

            delr = exptdelta(ranger, keys, spec, pos)


            rightChild = split + 1
            rightChildIsLeaf = (rightChild == ranger)
            #Threads.atomic_fence() # uncertain what this does and if it is necessary

            if !rightChildIsLeaf
                rightChild = branch_index(rightChild, spec)
            end


            boundingmax = max.(expand_integer(keys[rightChild].max), boundingmax)

            boundingmin = min.(expand_integer(keys[rightChild].min), boundingmin)
        else

            split = rangel - 1
            #rangel = Threads.atomic_cas!(store[split], K(0), ranger) 

            rangeleval = Atomix.@atomicreplace store[split] K(0) => ranger
            rangel = rangeleval[1]


            if rangel == 0 
                break
            end

            dell = exptdelta(rangel-1, keys, spec, pos)

            leftChild = split

            leftChildIsLeaf = (leftChild == rangel)
            
            #unclear if this is necessary
            #Threads.atomic_fence()
            if !leftChildIsLeaf
                leftChild = branch_index(leftChild, spec)
            end




            boundingmax = max.(expand_integer(keys[leftChild].max), boundingmax)

            boundingmin = min.(expand_integer(keys[leftChild].min), boundingmin)
        end


        q = delr < dell ? (ranger) : rangel 
         
        parentNode = branch_index(q, spec)
        #keys[parentNode].left = K(leftChild)
        pleftchild = leftChild
        #keys[parentNode] = GridKey{T,K}(keys[parentNode].index, keys[parentNode].morton_code, keys[parentNode].min, keys[parentNode].max, pleftchild, keys[parentNode].skip)
        keys[parentNode] = GridKey{T,K}(keys[parentNode].min, keys[parentNode].max, pleftchild, keys[parentNode].skip)
        

        if ranger == nL
            #keys[parentNode].skip = K(0)
            zerskip = K(0)
            keys[parentNode] = GridKey{T,K}(keys[parentNode].min, keys[parentNode].max, keys[parentNode].left, zerskip)
            #keys[parentNode] = GridKey{T,K}(boundingmin, boundingmax, keys[parentNode].left, zerskip)
        else
            r = ranger + K(1)
            if delr < exptdelta(r, keys, spec, pos)
                #keys[parentNode].skip = r
                keys[parentNode] = GridKey{T,K}(keys[parentNode].min, keys[parentNode].max, keys[parentNode].left, r)
                #keys[parentNode] = GridKey{T,K}(boundingmin, boundingmax, keys[parentNode].left, r)
            else
                #keys[parentNode].skip = branch_index(r, spec) 
                r = branch_index(r, spec) 
                keys[parentNode] = GridKey{T,K}(keys[parentNode].min, keys[parentNode].max, keys[parentNode].left, r)
                #keys[parentNode] = GridKey{T,K}(boundingmin, boundingmax, keys[parentNode].left, r)
            end
        end

        # keys[parentNode].min = boundingmin
        # keys[parentNode].max = boundingmax

        #keys[parentNode] = GridKey{T,K}(keys[parentNode].index, keys[parentNode].morton_code, boundingmin, boundingmax, keys[parentNode].left, keys[parentNode].skip)
        keys[parentNode] = GridKey{T,K}(compress_position(boundingmin, RoundDown), compress_position(boundingmax, RoundUp), keys[parentNode].left, keys[parentNode].skip)
        i = branch_index(q, spec)


        if i == branch_index(K(1), spec)
            return
        end
    end
end


function exptbounding_volume_hierarchy!(keys::Vector{GridKey{T,K}}, store, spec::SpheresBVHSpecs{T, K}, pos) where {T, K}

    @batch for i in 1:spec.leaves_count #in perfect parallel
        i = K(i)
        @inbounds exptbvh_interior!(keys, store, i, spec.leaves_count, spec.branches_count, spec, pos)
    end

    min = SVector{3, T}(0.0, 0.0, 0.0)
    max = SVector{3, T}(1.0, 1.0, 1.0)
    keys[branch_index(1, spec)] = GridKey{T,K}(min, max, keys[branch_index(1, spec)].left, keys[branch_index(1, spec)].skip)

end

###### Phase 3: traversal


Base.@propagate_inbounds function proximity_test!(neighbors::Vector{Tuple{K, K, T}},  query::IPointPrimitive{T,K}, subjects, spec::SpheresBVHSpecs{T, K}, squared_radius) where {T, K}
    #a2 = query.position .^ 2
    for each in eachindex(subjects)


        if query.index < subjects[each].index
            #dxyz2 = sum(a2 - 2 .* query.position .* subjects[each].position + (subjects[each].position .^ 2))
            dxyz2 = sum( (query.position - subjects[each].position) .^ 2 )

            #not demonstrably faster 3/15/2025
            #dxyz2 = sqeuclidean(positions[query_index].position, positions[leafsatoms].position) 

            # only push! new pairs that are close together 
            if dxyz2 < squared_radius

                

                d2 = sqrt( dxyz2 )
                #push!(neighbors, tuple(positions[query_index].index, positions[leafsatoms].index, dxyz, d2))
                push!(neighbors, tuple(query.index, subjects[each].index, d2))
            end

        end
    end

    return neighbors
end
Base.@propagate_inbounds function overlap_test(keys, currentKey, query_index, positions, spec::SpheresBVHSpecs{T, K}) where {T, K}
    return all(keys[currentKey].min .< positions[query_index].position .< keys[currentKey].max)

end

Base.@propagate_inbounds function exptoverlap_test(keys, currentKey, query_index, positions, spec::SpheresBVHSpecs{T, K}) where {T, K}
    #does this incure double getindex? could the cost be removed by tupling up the two friends?
    a = expand_integer(keys[currentKey].min)
    b = expand_integer(keys[currentKey].max)

    x = positions[query_index].position
    y = min.(max.(positions[query_index].position, a), b)
    return abs(sum(y .^ 2 ) - sum(x .^ 2)) < spec.neighbor_distance^2 #+ T(0.1)

    #return all(a .< positions[query_index].position .< b)

end

Base.@propagate_inbounds function neighbor_traverse(keys::Vector{GridKey{T,K}}, positions::Vector{IPointPrimitive{T,K}}, spec::SpheresBVHSpecs{T, K}) where {T, K}

    #neighbor_vec = parallel_neighbor_buffer(spec)
    
    threads = K(Threads.nthreads())

    #TODO isn't this structure extremely unfriendly to growing the individual elements at different times to different extents?
    #and what prevents death by thread racing??? certainly, nothing that i have consciously wrote
    # ----- not sure, having this as a vector of references seemed to diminish performance
    #neighbor_vec = [Vector{Tuple{K, K, SVector{3, T}, T}}(undef, 0) for i in 1:threads]
    neighbor_vec = parallel_neighbor_buffer(spec)
    #neighbor_vec = Vector{Vector{Tuple{K, K, T}}}(undef, Threads.nthreads())
    squared_radius = (spec.neighbor_distance) ^ 2


    @batch for chunk in 1:threads
         for query_index in K(chunk):threads:K(length(positions))
            currentKey = branch_index(1, spec)

            while currentKey != 0 # currentKey is the sentinel, end traversal of the given query
                # does query at all overlap with the volume of currentKey
                overlap = @inbounds overlap_test(keys, currentKey, query_index, positions, spec)
                # overlap = overlap_test(keys, currentKey, query_index, positions, spec)
                #overlap = sum(keys[currentKey].min .< positions[query_index].position .< keys[currentKey].max)
                # overlap = @btime overlap_test($keys, $currentKey, $query_index, $positions, $spec)
                # overlap = @btime newoverlap_test($keys, $currentKey, $query_index, $positions, $spec)
                # return

                if overlap

                    if keys[currentKey].left == 0 # currentKey is a leaf node

                        # i like the new method better, but i cannot for the life of me tell if one is better than the other. 
                        #when i repeat the same eval on the same data a million times, i see about a 10% improvement, but the behavior is def broken
                        # also, i like this new method better. Much prettier.
                        #@inbounds newproximity_test!(neighbor_vec[chunk], query_index, currentKey, positions, spec, squared_radius)
                        low = (currentKey-1) * spec.atomsperleaf + 1
                        high = currentKey * spec.atomsperleaf
                        slice = @view positions[low:high]
                        
                        @inbounds proximity_test!(neighbor_vec[chunk], positions[query_index], slice,  spec, squared_radius)

                        #proximity_test!(neighbor_vec[chunk], query_index, currentKey, positions, spec)
                        currentKey = keys[currentKey].skip
                    else #currentKey is a branch node, traverse to the left
                        currentKey = keys[currentKey].left
                    end
                else #query is not contained, can cut off traversal on the 'lefts' sequencef
                    currentKey = keys[currentKey].skip
                end

            end
        end
    end
    #return 
    #neighbors2 = neighbor_vec[1]
    for each in 2:1:threads
        append!(neighbor_vec[1], neighbor_vec[each] )
    end

    return neighbor_vec[1]

    # #single thread
    # # for each in eachindex(positions)
    # #     println(positions[each])
    # # end

    # for query_index in K(1):K(length(positions))#range(start=K(1), stop=K(spec.branches_count))
    #     currentKey = branch_index(1, spec)


    #     while currentKey != 0 # currentKey is the sentinel, end traversal of the given query
    #         # does query at all overlap with the volume of currentKey
    #         overlap = overlap_test(keys, currentKey, query_index, positions, spec)



    #         # if query_index==3
    #         #     println("Yees")
    #         #     println(currentKey)
    #         # end
    #         # if positions[query_index].index ==4 #&& positions[leafsatoms].index ==8
    #         #     println("query_index $query_index ")
    #         #     println("currentKey $currentKey ")
    #         #     #println("leafsatoms $leafsatoms ")
    #         #     println()
    #         # end

    #         if overlap == 3
    #             #println(keys[currentKey].left)
    #             if keys[currentKey].left == K(0) # currentKey is a leaf node
    #                 proximity_test!(neighbors, query_index, currentKey, positions, spec)
    #                 #proximity_test!(neighbors, query_index, keys[currentKey].index, positions, spec)
    #                 currentKey = keys[currentKey].skip
    #             else #currentKey is a branch node, traverse to the left
    #                 currentKey = keys[currentKey].left
    #             end
    #         else #query is not contained, can cut off traversal on the 'lefts' sequencef
    #             currentKey = keys[currentKey].skip
    #         end
    #     end
    # end

    return neighbors
end

function parallel_neighbor_buffer(spec::SpheresBVHSpecs{T, K}) where {T, K}
    #return [Vector{Tuple{K, K, SVector{3, T}, T}}(undef, 0) for i in 1:Threads.nthreads()]
    return [Vector{Tuple{K, K, T}}(undef, 0) for i in 1:Threads.nthreads()]
end

Base.@propagate_inbounds function expt_neighbor_traverse(keys::Vector{exptGridKey{T,K}}, positions::Vector{IPointPrimitive{T,K}}, spec::SpheresBVHSpecs{T, K}) where {T, K}

    #neighbor_vec = parallel_neighbor_buffer(spec)
    
    threads = K(Threads.nthreads())

    #TODO isn't this structure extremely unfriendly to growing the individual elements at different times to different extents?
    #and what prevents death by thread racing??? certainly, nothing that i have consciously wrote
    # ----- not sure, having this as a vector of references seemed to diminish performance
    #neighbor_vec = [Vector{Tuple{K, K, SVector{3, T}, T}}(undef, 0) for i in 1:threads]
    neighbor_vec = parallel_neighbor_buffer(spec)
    #neighbor_vec = Vector{Vector{Tuple{K, K, T}}}(undef, Threads.nthreads())
    squared_radius = (spec.neighbor_distance) ^ 2


    @batch for chunk in 1:threads
         for query_index in K(chunk):threads:K(length(positions))
            currentKey = branch_index(1, spec)

            while currentKey != 0 # currentKey is the sentinel, end traversal of the given query
                # does query at all overlap with the volume of currentKey
                overlap = @inbounds exptoverlap_test(keys, currentKey, query_index, positions, spec)
                overlapb = @inbounds overlap_test(keys, currentKey, query_index, positions, spec)
                # overlap = overlap_test(keys, currentKey, query_index, positions, spec)
                #overlap = sum(keys[currentKey].min .< positions[query_index].position .< keys[currentKey].max)
                # overlap = @btime overlap_test($keys, $currentKey, $query_index, $positions, $spec)
                # overlap = @btime newoverlap_test($keys, $currentKey, $query_index, $positions, $spec)
                # return
                #println(overlap==overlapb)
                if overlap
                    if keys[currentKey].left == 0 # currentKey is a leaf node

                        # i like the new method better, but i cannot for the life of me tell if one is better than the other. 
                        #when i repeat the same eval on the same data a million times, i see about a 10% improvement, but the behavior is def broken
                        # also, i like this new method better. Much prettier.
                        #@inbounds newproximity_test!(neighbor_vec[chunk], query_index, currentKey, positions, spec, squared_radius)
                        low = (currentKey-1) * spec.atomsperleaf + 1
                        high = currentKey * spec.atomsperleaf
                        slice = @view positions[low:high]
                        
                        @inbounds proximity_test!(neighbor_vec[chunk], positions[query_index], slice,  spec, squared_radius)

                        #proximity_test!(neighbor_vec[chunk], query_index, currentKey, positions, spec)
                        currentKey = keys[currentKey].skip
                    else #currentKey is a branch node, traverse to the left
                        currentKey = keys[currentKey].left
                    end
                else #query is not contained, can cut off traversal on the 'lefts' sequencef
                    currentKey = keys[currentKey].skip
                end

            end
        end
    end
    #return 
    #neighbors2 = neighbor_vec[1]
    for each in 2:1:threads
        append!(neighbor_vec[1], neighbor_vec[each] )
    end

    return neighbor_vec[1]

    #TODO better on data performance, worse on wall time, but also,, not?
    # somteims the profiler says this takes 20% of the time, and then other times it does nto!
    #return reduce(vcat, neighbor_vec)

    #return neighbors
end

###### Phase 4: put it all together
"""
    function build_traverse_bvh(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}
Build and then traverse a BVH, returning only a neighbor list
"""
function build_traverse_bvh(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}
    treeData = TreeData(position, spec)
    #neighbors = neighbor_traverse(treeData.tree, treeData.position, spec)
    #return neighbors
    return @inbounds neighbor_traverse(treeData[1], treeData[2], spec)
end

function exptbuild_traverse_bvh(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}
    treeData = exptTreeData(position, spec)
    #return expt_neighbor_traverse(treeData.tree, treeData.position, spec)
    return @inbounds expt_neighbor_traverse(treeData[5], treeData[2], spec)
end