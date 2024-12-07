export
    SpheresBVHSpecs,
    AxisAlignedBoundingBox,
    AABB,
    QuantizedAABB,
    GridKey,
    update_mortoncodes!,
    update_bvh!,
    build_bvh,
    traverse_bvh,
    batch_build_traverse,
    batched_batch_build
    #Atomic

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
    #this would be better, but it's easier to make something work with MVecs
    #maybe if the rest of this magically turns out to be realllly easy, i'll work here some more
    #min::Tuple{T, T, T}
    #max::Tuple{T, T, T}
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
    #parent_INode::T # this may be repurposed as the index of the skip connection
    #left::Ref{Union{GridKey, NaiveNode}}
    #right::Ref{Union{GridKey, NaiveNode}}
    left::T
    skip::T
end

##### Deprecated
# abstract type NaiveNode end
# mutable struct INode{T, K} <: NaiveNode
#     leaf_indices::Tuple{T, T}
#     #left::Ref{Union{GridKey, NaiveNode}}
#     #right::Ref{Union{GridKey, NaiveNode}}
#     left::T
#     skip::T
#     visits::T #to be marked as an @atomic, somehow
#     min::MVector{3, K}
#     max::MVector{3, K}
#     parent_INode::T
    
# end

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

function old_update_mortoncodes!(L, quantized_aabbs, morton_length, morton_type) 
    #TODO clean up this function to use L, aabbs, and spec?

    t3 = morton_type(3)
    t1 = morton_type(1)
    #L is an array of grid keys with an 'index' field which points to the 'nth index of a vector in the objectCollection struct' 
        # or a particular atom
    #n is the current nth bit of our morton code to change we wish to change, and it corresponds with every 3rd bit of our grid positions

    for each in eachindex(L)
        for n in morton_length 
            for dim in eachindex(quantized_aabbs[1].centroid)
                n_xyz = t3 * n + dim - t3 # we shift which position where on the morton code a bit should be written based on there being 3 spatial dimensions

                if (quantized_aabbs[each].centroid[dim]>>(n_xyz-t1)) & t1 != 0 #if a bit at index=n_xyz of a grid position is 1, then set the n_xyz'th bit of L[each] to one


                    L[each].morton_code ⊻= t1<<(n_xyz-t1)
                    

                else   
                    L[each].morton_code &= ~(t1<<(n_xyz-t1)) #set the n_xyz'th bit of L[each] to zero
                end

            end

        end


    end
end

##### this is the best version of the function.
# Curiously, it seems to be ~3x slower than the above version
# but at least this version produces the correct result

#some how this is 3x slower than the original version at half the operations and same allocation
"""
    function update_mortoncodes!(L, quantized_aabbs, morton_length, morton_type)

Take an array of GridKeys, L, an array of 3D integer coordinates, quantized_aabbs, and specification information,
to generate morton_codes for each GridKey.


"""
function update_mortoncodes!(L, quantized_aabbs, morton_length, morton_type) 

    inbit = morton_type(0)
    t3 = morton_type(3)
    t1 = morton_type(1)
    #L is an array of grid keys with an 'index' field which points to the 'nth index of a vector in the objectCollection struct' 
        # or a particular atom
    #n is the current nth bit of our morton code to change we wish to change, and it corresponds with every 3rd bit of our grid positions
    for each in eachindex(L)
        # set morton code to zero or else it will not work on second update call
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

function assign_mortoncodes(aabb_array, spec::SpheresBVHSpecs{T, K}, clct) where {T, K}
    morton_length = 0::Int64
    morton_string = " "::String
    #morton_type = K
    #idk how this should be better done to prevent runtime evaluation of a stupid if statement


    quantized_aabbs = [QuantizedAABB{K}(i, MVector{3, K}(0, 0, 0),  MVector{3, K}(0, 0, 0), MVector{3, K}(0, 0, 0)) for i in 1:spec.leaves_count]::Vector{QuantizedAABB{K}}

    update_gridkeys!(quantized_aabbs, aabb_array, spec, K, clct)

    L = [GridKey{K, T}(quantized_aabbs[i].index, 0, aabb_array[i].min, aabb_array[i].max, 0, 0) for i in 1:spec.leaves_count]

    update_mortoncodes!(L, quantized_aabbs, spec.morton_length, K)


    return L

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
    L = assign_mortoncodes(aabb_array, spec, clct)



    sort_mortoncodes!(L)




    
    return L

end



function update_bvh!(L, I, position, spec::SpheresBVHSpecs{T}, clct::Collector, aabb_array) where T
    update_aabb!(position, spec, aabb_array)
end


###### PHase 4: moving on to Prokopenko
# From Apetrei 2014 / Prokopenko and Lebrun-Grandie 2024, instead of countering common bits, were find the highest differing bit instead. 
## Though I am unconvinced this is the highest and not the lowest differing bit. Should be fine!


# according to Apetrei 2014, just returning the XOR is sufficient, finding the particular index of the relvant bit is unnecessary!

function delta_leaf(i, L::Vector{GridKey{K, T}}, spec::SpheresBVHSpecs{T,K}) where {T, K}

    # maintain numinternal nodes as length(L)-1 for now bc i dont think itll make a difference
    if  i >= spec.leaves_count || i < 1 #|| j > length(L) || j < 1 # i dont know if the same operation should be done to i, but idk how else to fix
        return typemax(K)
    end

    # XOR the two numbers to get a number with bits set where a and b differ
    # transcribed from ArborX directly, treeconstruction
    x = xor(L[i].morton_code, L[i+1].morton_code)


    
    return x + (x == K(0)) * (typemin(K) + (K(i) ⊻ K(i+1))) - K(1)


end
function delta_branch(i, L::Vector{GridKey{K, T}}, spec::SpheresBVHSpecs{T,K}) where {T, K}

    # maintain numinternal nodes as length(L)-1 for now bc i dont think itll make a difference
    if  i >= spec.leaves_count || i < 1 #|| j > length(L) || j < 1 # i dont know if the same operation should be done to i, but idk how else to fix
        return typemax(K)
    end

    # XOR the two numbers to get a number with bits set where a and b differ
    # transcribed from ArborX directly, treeconstruction
    x = xor(L[i].morton_code, L[i+1].morton_code)


    
    return x + (x == K(0)) * (typemin(K) + (K(i) ⊻ K(i+1))) - K(1)


end


        # if i == 2
        #     #println(L[rangel].morton_code, " icode")
        #     #println(L[rangel].morton_code, " rangelcode")
        #     #println(L[ranger].morton_code, " rangercode")
        #     println(rangel," ", ranger," ", " range l, r ")
        #     println(dell, " dell")
        #     println(delr, " delr")
        #     println()
        # end

function branch_index(a, spec::SpheresBVHSpecs{T, K}) where {T, K}
    return c = a + spec.leaves_count # this could do without the c = but i am not ready with testing yet
end

function stackless_reinterior!(store::Vector{Base.Threads.Atomic{Int64}}, i, nL, nI, keys, spec::SpheresBVHSpecs{T, K}) where {T, K}
    rangel = i 
    ranger = i
    dell = delta_leaf(rangel - 1, keys, spec)
    delr = delta_leaf(ranger, keys, spec)




    # stackless lBVH requires that the 'right most' elements of the tree point to an artificial node
    if i == nL # 1-based indexing requires this to be nL, but 0-based would be number of internal nodes
        keys[i].skip = 0
    else
        # Does Leaf[i+1] have a more similar morton code to Leaf[i] or Leaf[i+2]?
        # If Leaf[i+1] and Leaf[i] are more similar,
        # Then Leaf[i+1] and Leaf[i] are the right and left children of the same internal node.
        ir = i + 1
        #TODO is the <= necessary, or is < good enough?
        if delr < delta_leaf(ir, keys, spec) # are i and i+1 siblings, or are i+1 and i+2 siblings?
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

            delr = delta_branch(ranger, keys, spec)

            #here is wehre boundary computation is performed.
            # memory has to sync here for data safety

            # if i == branch_index(1, spec)
            #     println("leftchild in delr < dell $leftChild")
            # end

            # if leftChild == rangel
            #     println("left child is a branch in delr < dell")
            #     leftchild = branch_index(leftChild, spec)
            # end

            #leftChild = split
            # if leftChild == rangel-1
            #     println("yeah we made it here")
            #     leftChild = branch_index(leftChild, spec)
            # end

        else

            split = rangel - 1
            rangel = Threads.atomic_cas!(store[split], 0, ranger) 


            if rangel == 0 
                break
            end

            dell = delta_branch(rangel-1, keys, spec)

            leftChild = split
            if leftChild == rangel
                #println("arrived here")
                leftChild = branch_index(leftChild, spec)

            end

            # if i == branch_index(1, spec)
            #     println("leftchild in delr > dell $leftChild")
            # end
        end


        q = delr < dell ? (ranger) : rangel # commented out because the expanded version is easier to debug
         
        parentNode = branch_index(q, spec)


        if rangel == q
            keys[parentNode].left = leftChild
        else

            keys[parentNode].left = leftChild#branch_index(leftChild, spec )#leftChild#

        end
        

        if ranger == nL
            keys[parentNode].skip = 0
        else
            r = ranger + 1
            if delr < delta_branch(r, keys, spec)
                keys[parentNode].skip = r
            else
                keys[parentNode].skip = branch_index(r, spec) 
            end
        end

        
        
        i = branch_index(q, spec)


        if i == branch_index(1, spec)
            return
        end
    end
end

#by convention, if left or skip are negative, then they are referring to the index of the Inode, and positive is index of Leaf
function stacklessbottom_bvh(keys, spec::SpheresBVHSpecs{T, K}) where {T, K}

    store = [Base.Threads.Atomic{Int64}(0) for i in 1:spec.branches_count]



    #Threads.@threads 
    Threads.@threads :dynamic for i in spec.leaves_count:-1:1 #spec.leaves_count:-1:1#1:spec.leaves_count #in perfect parallel
        stackless_reinterior!(store, i, spec.leaves_count, spec.branches_count, keys, spec)


    end


end


# function sum_multi_good(a)
#     chunks = Iterators.partition(a, length(a) ÷ Threads.nthreads())
#     tasks = map(chunks) do chunk
#         Threads.@spawn sum_single(chunk)
#     end
#     chunk_sums = fetch.(tasks)
#     return sum_single(chunk_sums)
# end

###### Phase 6: traversal
# those times when you have to traverse down every path
# total sum of stack overflows / runaway Julia before I fixed this function and all its variants
# and caught recursion generating data conditions: 0
# to adequately prevent recursion, we would need a companion array for both the leaves and each of their connections
# as in, evaluate 'has this rope connection been used to update the traversal path already'
# and that's maybe too much work, so let's just run the above counter anyway.
function recursive_traversal(index, keys::Vector{GridKey{K, T}}, wasVisited::Vector{Bool}, spec::SpheresBVHSpecs{T, K}) where {T, K}
    traversal_count = 0
    left = keys[index].left
    skip = keys[index].skip
    wasVisited[index] = true
    # key about this structure: traverse TO the sentinel, then update values, then exit evaluation
    while true
        traversal_count += 1
        # no single traversal should run more times than contacting every branch and every leaf
        # but this won't stop runaway recursion, it will only stop a single call from creating infinite copies
        # however, there is currently the opportunity to make many recursion calls BEFORE exiting
        # overall, we will probably have infinitely growing recursion--especially if the tree has a loop
        if traversal_count > (spec.leaves_count + spec.branches_count)
            return
        end
        

        if index == 0 #||  # is the sentinel
            return
        else
            wasVisited[index] = true
            if (skip == 0) && (left == 0) # is pre-sentinel node, probably last leaf node
                index == skip

            elseif (skip == 0) && (left != 0)# is a right most node
                wasVisited[index] = true # TODO should be unnecessary now
                index = left
                skip = keys[left].skip
                left = keys[left].left              
                
            elseif (left == 0) && (skip != 0) #is a leaf node
                #do something
                wasVisited[index] = true
                index = skip
                left = keys[skip].left
                skip = keys[skip].skip
                
                
            elseif (skip != 0) && (left != 0) # is non-right branch, need to spawn a new traversal
                wasVisited[index] = true
                index = left
                recursive_traversal(skip, keys, wasVisited, spec)
                skip = keys[left].skip
                left = keys[left].left
                
                
            else
                println("unexpected traversal exit, datadump:")
                println("index $index, left $left, skip $skip")
                return 
            end
        end
    end
end



"""
    function is_traversable(keys::Vector{GridKey{K, T}}, spec::SpheresBVHSpecs{T, K}; ShowLonelyKeys=false) where {T, K})

Return a tuple of (Bool, "ReasonForResult"). Will print index of all unvisited keys if ShowLonelyKeys=true. Calls recursive_traversal which
calls itself at every internal branch that does not skip rope to the sentinel node.
"""
    
function is_traversable(keys::Vector{GridKey{K, T}}, spec::SpheresBVHSpecs{T, K}; ShowLonelyKeys=false) where {T, K}
    # should probably be an int to detect fro multiple visits

    # catch key pointing to itself, but not structural loops
    for each in eachindex(keys)
        if keys[each].left == each || keys[each].skip == each
            return false, "Self-referential"
        end
    end

    wasVisited = [false for each in eachindex(keys)]
    i = branch_index(1, spec)
    recursive_traversal(i, keys, wasVisited, spec)

    if ShowLonelyKeys
        for each in eachindex(wasVisited)
            if wasVisited[each] == false
                print("$each, ")
            end
        end
        println()
    end

    if sum(wasVisited) == length(keys)
        return true, "Fully traversable"
    else 
        return false, "Nontraversable"
    end

end
###### Phase 6.5: query traversal

###### Phase 7: put it all together
"""
    batch_build_traverse(runs::Int, position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T}) where {T, K}
Run a ```runs``` number of iterations of tree construction and recursive traversal on the same position data. Currently just collects
the number of good trees, self referential trees, or bad trees. My be expanded to show different tree variations if I get too bored.
"""

function batch_build_traverse(runs::Int, position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T}; printValues=true, printARun=false) where {T, K}
    goodTrees = 0
    selfishTrees = 0
    badTrees = 0

    for each in 1:runs
        keys = create_mortoncodes(position, spec, clct)::Vector{GridKey{K, T}} 
        I = [GridKey{K, T}(0, 0, MVector{3, K}(0.0, 0.0, 0.0), MVector{3, K}(0.0, 0.0, 0.0), 0, 0) for i in 1:spec.branches_count]
        append!(keys, I)
        stacklessbottom_bvh(keys, spec)
        result = is_traversable(keys, spec)
        if result[1]
            goodTrees += 1
        elseif result[2] == "Self-referential"
            selfishTrees += 1
        elseif result[2] == "Nontraversable"
            badTrees += 1
        end
    end
    if printValues
        println("goodTrees $goodTrees")
        println("selfishTrees $selfishTrees")
        println("badTrees $badTrees")
    else
        return goodTrees, selfishTrees, badTrees
    end

    if printARun 
        keys = create_mortoncodes(position, spec, clct)::Vector{GridKey{K, T}} 
        I = [GridKey{K, T}(0, 0, MVector{3, K}(0.0, 0.0, 0.0), MVector{3, K}(0.0, 0.0, 0.0), 0, 0) for i in 1:spec.branches_count]
        append!(keys, I)
        stacklessbottom_bvh(keys, spec)
        printtree(keys, spec)
    end
end

function batched_batch_build(setsOfRuns, runs::Int, position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T}) where {T, K}
    totalGood = 0
    totalSelfish = 0
    totalBad = 0
    for i in 1:setsOfRuns
        a = batch_build_traverse(runs, position, spec, clct; printValues=false, printARun=false)
        totalGood += a[1]
        totalSelfish += a[2]
        totalBad += a[3]
    end
    sum = totalGood + totalSelfish + totalBad
    println("% good ", totalGood / sum * 100)
    println("% selfish ", totalSelfish / sum * 100)
    println("% bad ", totalBad / sum * 100)
end

function printtree(keys::Vector{GridKey{K, T}}, spec::SpheresBVHSpecs{T, K}) where {T, K}
     println("Leaf connections")
    for i in 1:spec.leaves_count
        println(i, " ", keys[i].left, " ", keys[i].skip)
    end
    println()
    println("Branch connections")
    for i in branch_index(1, spec):1:branch_index(spec.branches_count, spec)
        println(i, " ", keys[i].left, " ", keys[i].skip)
    end
    println()
end
    


#TODO change this name?
function build_bvh(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T}) where {T, K}
    #L is an array of leaf nodes, 1 leaf per atom or 'primitive'
    keys = create_mortoncodes(position, spec, clct)::Vector{GridKey{K, T}} 
    ##println(L)
    #change type of int here to spec.morton_int, probably with Tuple(spec.morton_int[i, length(L)])

    #I = [tuple(i, length(L), Ref(L, i)[], Ref(L, i)[]) for i in 1:length(L)-1] # -1 from L because we want to have nodes = # atoms - 1 in this construction. Howard et al. chose a fixed 1024-1, but eh
    #I= [INode{K, T}(tuple(i, length(L)), 0, 0, 0,  MVector{3, T}(0, 0, 0), MVector{3, T}(0, 0, 0), 0) for i in 1:length(L)-1]::Vector{INode{K, T}} 
    I = [GridKey{K, T}(0, 0, MVector{3, K}(0.0, 0.0, 0.0), MVector{3, K}(0.0, 0.0, 0.0), 0, 0) for i in 1:spec.branches_count]
    append!(keys, I)

    #i = 3
    #j = 4
    ##println(bitstring(L[i].morton_code))
    ##println(bitstring(L[j].morton_code))
    #δ(L, i, j, spec)
    ##println(I)
    ##println(I)

    #bvh_solver!(L, I, spec)


    stacklessbottom_bvh(keys, spec)
    #println(is_traversable(keys, spec, ShowLonelyKeys=true))


    # counts1 = 0
    # counts4 = 0
    # counts5 = 0
    # countsDiff = 0
    # for i in 1:10000
    #     keys = create_mortoncodes(position, spec, clct)::Vector{GridKey{K, T}} 
    #     I = [GridKey{K, T}(0, 0, MVector{3, K}(0.0, 0.0, 0.0), MVector{3, K}(0.0, 0.0, 0.0), 0, 0) for i in 1:spec.branches_count]
    #     append!(keys, I)
    #     stacklessbottom_bvh(keys, spec)
    #     if keys[4].left == 4
    #         counts4 += 1
    #     elseif keys[4].left == 5
    #         counts5 += 1
    #     elseif keys[4].left == 1
    #         counts5 += 1
    #     else
    #         println("Weird keys[4]", keys[4].left)
    #     end
    # end
    # println("counts4 $counts4")
    # println("counts5 $counts5")
    # println("counts1 $counts1")
    # println("countsDiff $countsDiff")



    # for i in eachindex(L)
    #     println(del(i, L, spec))
    # end

    # for i in eachindex(position)
    #     println(position[i])
    # end
    printtree(keys, spec)

    
    
    #boundaries_wrapper(L, I, spec)
    #for i in eachindex(I)
        ##println("parentINode ",  I[i].parent_INode)#, ", ", I[i].max)
       # #println(I[i].min)
        ##println(I[i].max)
    #end
    ##println(I[1])
    ##println()
    #println(I[2].leaf_indices)
    #return neighborlist = println(traverse_bvh1(position, L, I, spec))
    #return traverse_bvh1(position, L, I, spec)
end