###### Bottom-up/parallel construction of linear bounding volume hierarchies with stackless traversal
# as described by Andrey Prokopenko and Damien Lebrun-Gradie
# https://arxiv.org/abs/2402.00665
# and as implemented by them in ArborX
# https://github.com/arborx/ArborX


# For tree generation, I tried to follow ArborX's implementation as closely as possible.
# Beyond that, I've made it up as I went along.




export
    SpheresBVHSpecs,
    MGridKey,
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
    mortoncodes!

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
    index::K
    morton_code::K
    min::SVector{3, T}
    max::SVector{3, T}
    left::K
    skip::K
end
mutable struct MGridKey{T, K}
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

# eachh point primitive has an 'index' to its original index place in the position (force and velocity) array(s)
##TODO expt if this can be immutable
mutable struct PointPrimitive{T, K}
    index::K
    morton_code::K #TODO update naming so this becomes one word, maybe even just morton?
    position::SVector{3, T}
end
"""
    struct TreeData{T, K}
        tree::Vector{GridKey{T, K}}
        position::Vector{IndexSafePosition{T, K}}
        quantizedposition::Vector{MVector{3, K}}
        store::Vector{Base.Threads.Atomic{Int64}}
    end
Generate a stackless bounding volume hierarchy stored in the `tree` field, while storing related data in tree construction for later reuse.
"""
struct TreeData{T, K}
    tree::Vector{MGridKey{T, K}}
    #position::Vector{IndexSafePosition{T, K}}
    position::Vector{PointPrimitive{T, K}}
    quantizedposition::Vector{MVector{3, K}}
    store::Vector{Base.Threads.Atomic{K}}
end
struct exptTreeData{T, K}
    tree::Vector{GridKey{T, K}}
    position::Vector{IndexSafePosition{T, K}}
    quantizedposition::Vector{MVector{3, K}}
    store::Vector{Base.Threads.Atomic{K}}

end

"""
    mortoncodes!(L, quantized_aabbs, spec::SpheresBVHSpecs{T, K}) where {T, K} )

Take an array of GridKeys, L, an array of 3D integer coordinates, quantized aabbs, and specification information,
to generate morton codes for each GridKey.

For a discussion into how this function works, see 11 January, 2025 in devdiary, approx. L1245.
"""
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
        L[each] = GridKey{T,K}(L[each].index, input, L[each].min, L[each].max, L[each].left, L[each].skip)
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
        L[each].morton_code = K(0) 
        result = K(0)
        for dim in eachindex(quantized[1])
            yin = quantized[each][dim] & magic_values[dim]
            yang = quantized[each][dim] | magic_not_values[dim]
            yinyang = yin & yang

            L[each].morton_code |= yinyang

        end
    end

end
function sort_mortoncodes!(L, spec::SpheresBVHSpecs{T, K}) where {T, K}
    #make a specialized radix sort to replace base sort // space for GPU backends to put forward their own sort
    sort!(L, by = x -> x.morton_code, rev=false) # sorts lexicographically both the binary and the integer
    #sort!(L, by=x -> count(c -> c == '1', bitstring(x.morton_code)))
    #sort!(L, by = x -> x.morton_code), alg=RadixSort #wont run, RadixSort does not have iterate defined
end



#TODO how can this work with stativ vectors instead, not for speed, just for Static vector unity
function quantized_positions!(quantized::Vector{MVector{3, K}}, pos::Vector{PointPrimitive{T,K}}, spec::SpheresBVHSpecs{T, K}) where {T, K}

    for dim in eachindex(quantized[1])
        sort!(pos, by=x->x.position[dim])
        for each in eachindex(pos)
            quantized[pos[each].index][dim] = K(each)
        end
    end
end

#TODO make MGridKey function?
function cluster_primitives(L::Vector{PointPrimitive{T, K}}, spec::SpheresBVHSpecs{T, K}) where {T, K}
    #leaves = Vector{MGridKey{T, K}}(undef, spec.leaves_count)
    leaves = [MGridKey{T,K}(0, 0, SVector{3, T}(0.0, 0.0, 0.0), SVector{3, T}(0.0, 0.0, 0.0), 0, 0) for i in 1:spec.leaves_count]
    for each in eachindex(leaves)#firstindex(L):spec.atomsperleaf:lastindex(L)
        #firs
        #println("each $each")
        #TODO fix this stuff
        #rabbit = MGridKey{T,K}(0, 0, SVector{3, T}(0.0, 0.0, 0.0), SVector{3, T}(0.0, 0.0, 0.0), 0, 0)
        bird = (each-1) * spec.atomsperleaf + 1
        println("spec leaves ", spec.leaves_count)
        println("bird $bird")
        println("length L ", length(L))
        leaves[each].min = L[bird].position .- spec.neighbor_distance
        leaves[each].max = L[bird].position .+ spec.neighbor_distance

        for i in 1:1:spec.atomsperleaf#spec.atomsperleaf-1:-1:1#0:1:spec.atomsperleaf - 1#+1 #this addition increases passes by 1, but i have no idea why
            #println("i $i ")
            a = (each-1)*spec.atomsperleaf  + i 

            #println("a ", a)
            leaves[each] = MGridKey{T,K}(L[each].index, 
                                    L[each].morton_code, 
                                    # min.(L[each].min, L[a].min), 
                                    # max.(L[each].max, L[a].max),
                                    min.(leaves[each].min, L[a].position .- spec.neighbor_distance),  # necessary to add neighbor distance when 1 atom per leaf
                                    max.(leaves[each].max, L[a].position .+ spec.neighbor_distance),
                                    0,
                                    0
            )
        end

    end

    return leaves
end

#TODO should this be capitalized??? im waffling
function exptTreeData(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}

    pos = [IndexSafePosition{T, K}(i, SVector{3, T}(position[i])) for i in eachindex(position)]
    #println(isbits(pos))

    quantized_xyz = [MVector{3, K}(0, 0, 0) for i in eachindex(position)]

    quantized_positions!(quantized_xyz, pos, spec)

    # newL = @time [GridKey{T, K}(i, 0, 
    #         SVector{3, T}(position[i] .- spec.bounding_distance), SVector{3, T}(position[i] .+ spec.bounding_distance),
    #         0, 0) for i in 1:spec.leaves_count
    # ]
    # println(length(newL))
    # println(spec.leaves_count)
    #newL = SizedVector{spec.leaves_count, GridKey{T,K}}(newL)
    #TODO really really heavy alloc pressure here
    L = [GridKey{T, K}(i, 0, 
            SVector{3, T}(position[i] .- spec.bounding_distance), SVector{3, T}(position[i] .+ spec.bounding_distance),
            0, 0) for i in 1:spec.leaves_count
    ]
    #println(isbits(L))
    #println(isbits(newL))
    exptmortoncodes!(L, quantized_xyz, spec)
    #@btime mortoncodes!($L, $quantized_xyz, $spec)
    #@btime mortoncodes_imm!($newL, $quantized_xyz, $spec)
    sort_mortoncodes!(L, spec)
    #@btime sort_mortoncodes!($L, $spec)
    #@btime sort_mortoncodes!($newL, $spec)
    store = [Base.Threads.Atomic{K}(0) for i in 1:spec.branches_count]



    I = [GridKey{T, K}(0, 0, SVector{3, T}(0.0, 0.0, 0.0), SVector{3, K}(0.0, 0.0, 0.0), 0, 0) for i in 1:spec.branches_count]
    
    append!(L, I)
    immuteL = [GridKey{T, K}(L[each].index, L[each].morton_code, L[each].min, L[each].max, L[each].left, L[each].skip) for each in eachindex(L)]

    exptbounding_volume_hierarchy!(L, store, spec)



    return exptTreeData{T, K}(L, pos, quantized_xyz, store)
end
function TreeData(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}


    # hackspec = SpheresBVHSpecs(; bounding_distance=spec.bounding_distance, 
    #                                 neighbor_distance=spec.neighbor_distance, 
    #                                 leaves_count=length(position) / spec.atomsperleaf,
    #                                 floattype=Float32,
    #                                 atomsperleaf = spec.atomsperleaf 
    # )

    #pos = [IndexSafePosition{T, K}(i, SVector{3, T}(position[i])) for i in eachindex(position)]
    #println(isbits(pos))
    L = [PointPrimitive{T,K}(i, 0, position[i]) for i in 1:spec.atom_count]

    quantized_xyz = [MVector{3, K}(0, 0, 0) for i in eachindex(position)]

    quantized_positions!(quantized_xyz, L, spec)

    # newL = @time [GridKey{T, K}(i, 0, 
    #         SVector{3, T}(position[i] .- spec.bounding_distance), SVector{3, T}(position[i] .+ spec.bounding_distance),
    #         0, 0) for i in 1:spec.leaves_count
    # ]
    # println(length(newL))
    # println(spec.leaves_count)
    #newL = SizedVector{spec.leaves_count, GridKey{T,K}}(newL)
    #TODO really really heavy alloc pressure here

    # L = [MGridKey{T, K}(i, 0, 
    # SVector{3, T}(position[i] .- spec.bounding_distance), SVector{3, T}(position[i] .+ spec.bounding_distance),
    # 0, 0) for i in 1:spec.leaves_count
    # ]
#TODO is treshold shift necessary here?
    # L = [MGridKey{T, K}(i, 0, 
    #         SVector{3, T}(position[i] .- spec.neighbor_distance), SVector{3, T}(position[i] .+ spec.neighbor_distance),
    #         0, 0) for i in 1:spec.leaves_count
    # ]
    



    #println(isbits(L))
    #println(isbits(newL))
    mortoncodes!(L, quantized_xyz, spec)
    #@btime mortoncodes!($L, $quantized_xyz, $spec)
    #@btime mortoncodes_imm!($newL, $quantized_xyz, $spec)
    sort_mortoncodes!(L, spec)
    #@btime sort_mortoncodes!($L, $spec)
    #@btime sort_mortoncodes!($newL, $spec)

    #store = [Base.Threads.Atomic{K}(0) for i in 1:spec.branches_count]
    store = [Base.Threads.Atomic{K}(0) for i in 1:spec.branches_count]
    #lenL = length(L)
    #hackL = Vector{MGridKey{T, K}}(undef, 0)
    hackL = cluster_primitives(L, spec)
    # for each in firstindex(L):spec.atomsperleaf:lastindex(L)
    #     #firs
    #     #println("each $each")
    #     #TODO fix this stuff
    #     rabbit = MGridKey{T,K}(0, 0, SVector{3, T}(0.0, 0.0, 0.0), SVector{3, T}(0.0, 0.0, 0.0), 0, 0)
    #     for i in 0:1:spec.atomsperleaf-1
    #         #println("i $i ")
    #         a = each+i
    #         if a > length(L)
    #             break
    #         end
    #         #println("a ", a)
    #         rabbit = MGridKey{T,K}(L[each].index, 
    #                                 L[each].morton_code, 
    #                                 # min.(L[each].min, L[a].min), 
    #                                 # max.(L[each].max, L[a].max),
    #                                 min.(L[each].position .- spec.neighbor_distance, L[a].position .- spec.neighbor_distance),  # necessary to add neighbor distance when 1 atom per leaf
    #                                 max.(L[each].position .+ spec.neighbor_distance, L[a].position .+ spec.neighbor_distance),
    #                                 0,
    #                                 0
    #         )
    #     end
    #     push!(hackL, rabbit)
    # end


    # I = [MGridKey{T, K}(0, 0, SVector{3, T}(0.0, 0.0, 0.0), SVector{3, K}(0.0, 0.0, 0.0), 0, 0) for i in 1:spec.branches_count]
    
    # append!(L, I)

    # bounding_volume_hierarchy!(L, store, spec)

    # return TreeData{T, K}(L, pos, quantized_xyz, store)


    hackI = [MGridKey{T, K}(0, 0, SVector{3, T}(0.0, 0.0, 0.0), SVector{3, K}(0.0, 0.0, 0.0), 0, 0) for i in 1:spec.branches_count]
    
    append!(hackL, hackI)

    bounding_volume_hierarchy!(hackL, store, spec)

    return TreeData{T, K}(hackL, L, quantized_xyz, store)
end



function TreeData!(treeData::TreeData{T, K}, position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}

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

    bounding_volume_hierarchy!(treeData.tree, treeData.store, spec)
    
end

###### Phase 2: tree construction


function delta(i, L, spec::SpheresBVHSpecs{T,K}) where {T, K}

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
    return K(a + spec.leaves_count)
end


function bvh_interior!(keys::Vector{MGridKey{T, K}}, store::Vector{Base.Threads.Atomic{K}}, i, nL, nI, spec::SpheresBVHSpecs{T, K}) where {T, K}
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
            ranger = Threads.atomic_cas!(store[split], K(0), rangel)



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

            #boundingmax = (keys[rightChild].max .< boundingmax) .* keys[rightChild].max .+ (boundingmax .< keys[rightChild].max) .* boundingmax
            boundingmax = max.(keys[rightChild].max, boundingmax)
            #TODO i think having this also is a correct read of arborx?
            boundingmin = min.(keys[rightChild].min, boundingmin)
        else

            split = rangel - 1
            rangel = Threads.atomic_cas!(store[split], K(0), ranger) 


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


            #boundingmin = (keys[leftChild].min .< boundingmin) .* keys[leftChild].min .+ (boundingmin .< keys[leftChild].min) .* boundingmin
            
            boundingmax = max.(keys[leftChild].max, boundingmax)
            #TODO i think having this also is a correct read of arborx?
            boundingmin = min.(keys[leftChild].min, boundingmin)

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


function bounding_volume_hierarchy!(keys, store, spec::SpheresBVHSpecs{T, K}) where {T, K}

    Threads.@threads for i in 1:spec.leaves_count #in perfect parallel
        i = K(i)
        bvh_interior!(keys, store, i, spec.leaves_count, spec.branches_count, spec)
    end

    keys[branch_index(1, spec)].min = SVector{3, T}(0.0, 0.0, 0.0)
    keys[branch_index(1, spec)].max = SVector{3, T}(1.0, 1.0, 1.0)

end

function exptbvh_interior!(keys, store::Vector{Base.Threads.Atomic{K}}, i, nL, nI, spec::SpheresBVHSpecs{T, K}) where {T, K}
    rangel = i 
    ranger = i
    dell = delta(rangel - 1, keys, spec)
    delr = delta(ranger, keys, spec)
    boundingmin = keys[i].min
    boundingmax = keys[i].max





    # stackless lBVH requires that the 'right most' elements of the tree point to an artificial node
    if i == nL # 1-based indexing requires this to be nL, but 0-based would be number of internal nodes
        #keys[i].skip = K(0)
        zerskip = K(0)
        keys[i] = GridKey{T,K}(keys[i].index, keys[i].morton_code, keys[i].min, keys[i].max, keys[i].left, zerskip)
    else
        # Does Leaf[i+1] have a more similar morton code to Leaf[i] or Leaf[i+2]?
        # If Leaf[i+1] and Leaf[i] are more similar,
        # Then Leaf[i+1] and Leaf[i] are the right and left children of the same internal node.
        ir = i + K(1)

        if delr < delta(ir, keys, spec) # are i and i+1 siblings, or are i+1 and i+2 siblings?
            #keys[i].skip = ir
            zerskip = K(ir)
            keys[i] = GridKey{T,K}(keys[i].index, keys[i].morton_code, keys[i].min, keys[i].max, keys[i].left, zerskip)
        else
            
            #keys[i].skip = branch_index(K(ir), spec)
            zerskip = branch_index(K(ir), spec)
            keys[i] = GridKey{T,K}(keys[i].index, keys[i].morton_code, keys[i].min, keys[i].max, keys[i].left, zerskip)
            
        end
    end

    while true
        isLeftChild = delr < dell
        if isLeftChild
            leftChild = i

            split = ranger # split position between the range of keys covered by any given INode
            ranger = Threads.atomic_cas!(store[split], K(0), rangel)



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
            rangel = Threads.atomic_cas!(store[split], K(0), ranger) 


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
        #keys[parentNode].left = K(leftChild)
        pleftchild = leftChild
        keys[parentNode] = GridKey{T,K}(keys[parentNode].index, keys[parentNode].morton_code, keys[parentNode].min, keys[parentNode].max, pleftchild, keys[parentNode].skip)

        

        if ranger == nL
            #keys[parentNode].skip = K(0)
            zerskip = K(0)
            keys[parentNode] = GridKey{T,K}(keys[parentNode].index, keys[parentNode].morton_code, keys[parentNode].min, keys[parentNode].max, keys[parentNode].left, zerskip)
        else
            r = ranger + K(1)
            if delr < delta(r, keys, spec)
                #keys[parentNode].skip = r
                keys[parentNode] = GridKey{T,K}(keys[parentNode].index, keys[parentNode].morton_code, keys[parentNode].min, keys[parentNode].max, keys[parentNode].left, r)
            else
                #keys[parentNode].skip = branch_index(r, spec) 
                r = branch_index(r, spec) 
                keys[parentNode] = GridKey{T,K}(keys[parentNode].index, keys[parentNode].morton_code, keys[parentNode].min, keys[parentNode].max, keys[parentNode].left, r)
            end
        end

        # keys[parentNode].min = boundingmin
        # keys[parentNode].max = boundingmax

        keys[parentNode] = GridKey{T,K}(keys[parentNode].index, keys[parentNode].morton_code, boundingmin, boundingmax, keys[parentNode].left, keys[parentNode].skip)
        
        i = branch_index(q, spec)


        if i == branch_index(K(1), spec)
            return
        end
    end
end


function exptbounding_volume_hierarchy!(keys::Vector{GridKey{T,K}}, store, spec::SpheresBVHSpecs{T, K}) where {T, K}

    Threads.@threads for i in 1:spec.leaves_count #in perfect parallel
        i = K(i)
        exptbvh_interior!(keys, store, i, spec.leaves_count, spec.branches_count, spec)
    end

    min = SVector{3, T}(0.0, 0.0, 0.0)
    max = SVector{3, T}(1.0, 1.0, 1.0)
    keys[branch_index(1, spec)] = GridKey{T,K}(keys[branch_index(1, spec)].index, keys[branch_index(1, spec)].morton_code, min, max, keys[branch_index(1, spec)].left, keys[branch_index(1, spec)].skip)

end

###### Phase 3: traversal
function proximity_test!(neighbors::Vector{Tuple{K, K, T}},  query_index::K, currentKey::K, positions::Vector{PointPrimitive{T,K}}, spec::SpheresBVHSpecs{T, K}) where {T, K}
    
    # eliminate redundant pairings when qi == cK 
    #AND when [qi = a, cK = b] with [qi = b, cK = a] appear in the same pairs list array
    # if !(query_index < currentKey) 
    #     return
    # else
    #     #TODO does dot syntax screw with thigns hwer?

    #         dxyz = (positions[query_index] - positions[currentKey])
    #         d2 = sqrt( sum(dxyz .^ 2))

    #         # only push! new pairs that are close together 
    #         if d2 < spec.neighbor_distance
    #             push!(neighbors, tuple(query_index, currentKey, dxyz, d2))

    #         end
    # end
    # if query_index==3
    #     println("Yees")
    #     println(currentKey)
    # end
    for i in 1:1:spec.atomsperleaf #- 1
        #localquery = query_index + i #-1 or smth. #TODO should be query index or current key? what part did i bung uyp here??
        #println(i)
        
        leafsatoms = (currentKey-1) * spec.atomsperleaf  + i 
        # if positions[query_index].index ==6 #&& positions[leafsatoms].index ==8
        #     println("query_index $query_index ")
        #     println("currentKey $currentKey ")
        #     println("leafsatoms $leafsatoms ")
        #     println()
        # end
        # if positions[query_index].index ==3 || positions[leafsatoms].index ==5
        #     println("query_index $query_index ")
        #     println("currentKey $currentKey ")
        #     println("leafsatoms $leafsatoms ")
        #     println()
        # end
        # if positions[query_index].index == 94  & positions[leafsatoms].index > 98
        #     println("query_index $query_index ")
        #     println("currentKey $currentKey ")
        #     println("leafsatoms $leafsatoms ")
        #     #println("leafatoms index ", positions[leafsatoms].index)
        #     println()
        # end

        #if localquery < currentKey 
        if positions[query_index].index < positions[leafsatoms].index #TODO as an if else

            #TODO does dot syntax screw with thigns hwer?

            #dxyz = (positions[localquery] - positions[currentKey])
            dxyz = (positions[query_index].position - positions[leafsatoms].position)
            d2 = sqrt( sum(dxyz .^ 2))

            # only push! new pairs that are close together 
            if d2 < spec.neighbor_distance

 
                #push!(neighbors, tuple(localquery, currentKey, dxyz, d2))
                #push!(neighbors, tuple(positions[query_index].index, positions[leafsatoms].index, dxyz, d2))
                push!(neighbors, tuple(positions[query_index].index, positions[leafsatoms].index, d2))
            end

        end
    end

    return neighbors
end
@inline function overlap_test(keys, currentKey, query_index, positions, spec::SpheresBVHSpecs{T, K}) where {T, K}
    return sum(keys[currentKey].min .< positions[query_index].position .< keys[currentKey].max)
end

@inline function newoverlap_test(keys::Vector{GridKey{T,K}}, currentKey::K, query_index::K, positions, spec::SpheresBVHSpecs{T, K}) where {T, K}

    return prod(keys[currentKey].min .< positions[query_index] .< keys[currentKey].max) 
end

function neighbor_traverse(keys::Vector{MGridKey{T,K}}, positions::Vector{PointPrimitive{T,K}}, spec::SpheresBVHSpecs{T, K}) where {T, K}

    # threads = K(Threads.nthreads())

    # #TODO isn't this structure extremely unfriendly to growing the individual elements at different times to different extents?
    # #and what prevents death by thread racing??? certainly, nothing that i have consciously wrote
    # # ----- not sure, having this as a vector of references seemed to diminish performance
    # #neighbor_vec = [Vector{Tuple{K, K, SVector{3, T}, T}}(undef, 0) for i in 1:threads]
    # neighbor_vec = [Vector{Tuple{K, K, T}}(undef, 0) for i in 1:threads]


    # Threads.@threads for chunk in 1:threads
    #     for query_index in K(chunk):threads:K(length(positions))#range(start=K(1), stop=K(spec.branches_count))
    #         currentKey = branch_index(1, spec)

    #         while currentKey != 0 # currentKey is the sentinel, end traversal of the given query
    #             # does query at all overlap with the volume of currentKey
    #             overlap = overlap_test(keys, currentKey, query_index, positions, spec)





    #             if overlap > 0
    #                 if keys[currentKey].left == 0 # currentKey is a leaf node
    #                     proximity_test!(neighbor_vec[chunk], query_index, currentKey, positions, spec)
    #                     currentKey = keys[currentKey].skip
    #                 else #currentKey is a branch node, traverse to the left
    #                     currentKey = keys[currentKey].left
    #                 end
    #             else #query is not contained, can cut off traversal on the 'lefts' sequencef
    #                 currentKey = keys[currentKey].skip
    #             end

    #         end
    #     end
    # end
    # neighbors = neighbor_vec[1]
    # for each in 2:1:threads
    #     append!(neighbors, neighbor_vec[each] )
    # end

    #neighbors = [Tuple{K, K, SVector{3, T}, T}(undef)]


    #neighbors = Vector{Tuple{K, K, SVector{3, T}, T}}(undef, 0)
    neighbors = Vector{Tuple{K, K, T}}(undef, 0)
    # for each in eachindex(keys)
    #     println(keys[each])
    # end

    #single thread
    # for each in eachindex(positions)
    #     println(positions[each])
    # end

    for query_index in K(1):K(length(positions))#range(start=K(1), stop=K(spec.branches_count))
        currentKey = branch_index(1, spec)


        while currentKey != 0 # currentKey is the sentinel, end traversal of the given query
            # does query at all overlap with the volume of currentKey
            overlap = overlap_test(keys, currentKey, query_index, positions, spec)



            # if query_index==3
            #     println("Yees")
            #     println(currentKey)
            # end
            # if positions[query_index].index ==4 #&& positions[leafsatoms].index ==8
            #     println("query_index $query_index ")
            #     println("currentKey $currentKey ")
            #     #println("leafsatoms $leafsatoms ")
            #     println()
            # end

            if overlap > 0
                #println(keys[currentKey].left)
                if keys[currentKey].left == K(0) # currentKey is a leaf node
                    proximity_test!(neighbors, query_index, currentKey, positions, spec)
                    #proximity_test!(neighbors, query_index, keys[currentKey].index, positions, spec)
                    currentKey = keys[currentKey].skip
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

function parallel_neighbor_buffer(spec::SpheresBVHSpecs{T, K}) where {T, K}
    return [Vector{Tuple{K, K, SVector{3, T}, T}}(undef, 0) for i in 1:Threads.nthreads()]
end

function expt_neighbor_traverse(keys::Vector{GridKey{T,K}}, positions::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}

    #threads = K(Threads.nthreads())

    #TODO isn't this structure extremely unfriendly to growing the individual elements at different times to different extents?
    # ----- not sure, having this as a vector of references seemed to diminish performance
    #neighbor_vec = [Vector{Tuple{K, K, SVector{3, T}, T}}(undef, 0) for i in 1:threads]
    neighbor_vec = parallel_neighbor_buffer(spec)
    #neighbor = Vector{Tuple{K, K, SVector{3, T}, T}}(undef, 0)


    # clamp traversal to 1 index before the last leaf because if every other leaf has been considered, 
    # then the last leaf does not need to be reconsidered for a neighbor pair
#     currentKey = 5
#     query_index = 2

#     overlap = @btime overlap_test($keys, $currentKey, $query_index, $positions, $spec)
#    # overlapb = @btime newoverlap_test($keys, $currentKey, $query_index, $positions, $spec)
#     println()
#     return
    

Threads.@threads for chunk in 1:Threads.nthreads()
        for query_index in chunk:Threads.nthreads():spec.branches_count#range(start=K(1), stop=K(spec.branches_count))
            currentKey = branch_index(1, spec)
            query_index = K(query_index)

            while currentKey != 0 # currentKey is the sentinel, end traversal of the given query
                # does query at all overlap with the volume of currentKey
                #overlap = newoverlap_test(keys, currentKey, query_index, positions, spec)
                overlap = overlap_test(keys, currentKey, query_index, positions, spec)





                if overlap > 0
                    if keys[currentKey].left == 0 # currentKey is a leaf node
                        proximity_test!(neighbor_vec[chunk], query_index, currentKey, positions, spec)
                        #proximity_test!(neighbor, query_index, currentKey, positions, spec)
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
    if Threads.nthreads() > 1
            
        neighbors = neighbor_vec[1]
        for each in 2:1:Threads.nthreads()
            append!(neighbors, neighbor_vec[each] )
        end

        return neighbors
    else
        return neighbor_vec[1]
    end
    #return neighbor
end

###### Phase 4: put it all together
"""
    function build_traverse_bvh(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}
Build and then traverse a BVH, returning only a neighbor list
"""
function build_traverse_bvh(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}
    #     hackspec = SpheresBVHSpecs(;  
    #                                     neighbor_distance=spec.neighbor_distance, 
    #                                     atom_count=length(position),
    #                                     floattype=Float32,
    #                                     atomsperleaf = spec.atomsperleaf 
    # )
    treeData = TreeData(position, spec)
    neighbors = neighbor_traverse(treeData.tree, treeData.position, spec)
    return neighbors
end

function exptbuild_traverse_bvh(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}
    treeData = exptTreeData(position, spec)
    return expt_neighbor_traverse(treeData.tree, position, spec)
end