module NaiveKA

using NaiveDynamics
using StaticArrays
using Adapt
using KernelAbstractions
using Atomix
using BenchmarkTools


###### Bottom-up/parallel construction of linear bounding volume hierarchies with stackless traversal
# as described by Andrey Prokopenko and Damien Lebrun-Gradie
# https://arxiv.org/abs/2402.00665
# and as implemented by them in ArborX
# https://github.com/arborx/ArborX


# For tree generation, I tried to follow ArborX's implementation as closely as possible.
# Beyond that, I've made it up as I went along.

# GridKey compression into bin coordinates was implemented as described in Howard et al. 2019
#https://doi.org/10.1016/j.commatsci.2019.04.004




# This GPU implementation was prepared thanks yo the documentation of KernelAbstractions.jl 
# and Tim Besard's talk: https://www.youtube.com/watch?v=QvlBmh1t9I4


##### Phase 1: morton code construction and updating

#TODO does this need to be mutable? Can we improve perf by restructuring to be immutable
# can index of the original atom be stored in a companion array? would this matter?
# struct GridKey{T, K}
#     #index::K
#     #morton_code::K
#     min::SVector{3, T}
#     max::SVector{3, T}
#     left::K
#     skip::K
# end
# struct exptGridKey{T, K}
#     min::Int32 #in a float64 scenario, Int64 would not capture anywhere near the ratio of significance that Int32 can capture from 3 float32s
#     max::Int32
#     left::K
#     skip::K
# end

# function Adapt.adapt_structure(to, spec::SpheresBVHSpecs)
#     neighbor_distance = adapt_structure(to, spec.neighbor_distance)
#     atom_count = adapt_structure(to, spec.atom_count)
#     leaves_count = adapt_structure(to, spec.leaves_count)
#     branches_count = adapt_structure(to, spec.branches_count)
#     atomsperleaf = adapt_structure(to, spec.atomsperleaf)
#     SpheresBVHSpecs(neighbor_distance, atom_count, leaves_count, branches_count, atomsperleaf)
# end

# function Adapt.adapt_structure(to, key::exptGridKey)
#     min = adapt_structure(to, key.min)
#     max = adapt_structure(to, key.max)
#     left = adapt_structure(to, key.left)
#     skip = adapt_structure(to, key.skip)
#     GridKey(min, max, left, skip)
# end

# struct PointPrimitive{T, K}
#     index::K
#     morton_code::K
#     position::SVector{3, T}
# end

Adapt.@adapt_structure NaiveDynamics.SpheresBVHSpecs
Adapt.@adapt_structure NaiveDynamics.GridKey

const binwidth = Float32(1/1023)
const KA = KernelAbstractions

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



function delta(i, pos, spec::SpheresBVHSpecs{T,K}) where {T, K}

    # maintain numinternal nodes as length(L)-1 for now bc i dont think itll make a difference
    if  i >= spec.leaves_count || i < 1 #|| j > length(L) || j < 1 # i dont know if the same operation should be done to i, but idk how else to fix
        return typemax(K)
    end

# eachh point primitive has an 'index' to its original index place in the position (force and velocity) array(s)
##TODO expt if this can be immutable
# mutable struct PointPrimitive{T, K}
#     index::K
#     morton_code::K #TODO update naming so this becomes one word, maybe even just morton?
#     position::SVector{3, T}
# end
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


function bvh_interior!(keys, store, i, nL, nI, spec::SpheresBVHSpecs{T, K}, pos) where {T, K}

end


KA.@kernel cpu = false function NaiveDynamics.gpubounding_volume_hierarchy!(keys::AbstractVector{GridKey{T,K}}, store::AbstractVector{K}, @Const(spec::SpheresBVHSpecs{T, K}), @Const(pos)) where {T, K}
    
    i = KA.@index(Global, Linear)
    
    #KA.@print(KA.@groupsize)
    
    rangel = i
    
    ranger = i
    dell = delta(rangel - 1, pos, spec) 
    delr = delta(ranger, pos, spec) 
    #KA.@print("happy")
    
    boundingmin = keys[i].min
    
    boundingmax = keys[i].max
    
    


    #KA.@print("happier")
    #KA.@print(i)
    #KA.@print(spec.leaves_count)
    # stackless lBVH requires that the 'right most' elements of the tree point to an artificial node
    if i == spec.leaves_count # 1-based indexing requires this to be nL, but 0-based would be number of internal nodes
        #keys[i].skip = K(0)
        
        zerskip = K(0)
        
        #keys[i] = GridKey{T,K}(keys[i].index, keys[i].morton_code, keys[i].min, keys[i].max, keys[i].left, zerskip)
        keys[i] = GridKey{T,K}(keys[i].min, keys[i].max, keys[i].left, zerskip)
        
    else
        # Does Leaf[i+1] have a more similar morton code to Leaf[i] or Leaf[i+2]?
        # If Leaf[i+1] and Leaf[i] are more similar,
        # Then Leaf[i+1] and Leaf[i] are the right and left children of the same internal node.
        ir = i + K(1)
        #KA.@print("happiest")
        if delr < delta(ir, pos, spec) # are i and i+1 siblings, or are i+1 and i+2 siblings?
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
            
            delr = delta(ranger, pos, spec)


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

            dell = delta(rangel-1, pos, spec)

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
        #KA.@print("happiest")
        

        q = delr < dell ? (ranger) : rangel 
            
        parentNode = branch_index(q, spec)
        #keys[parentNode].left = K(leftChild)
        pleftchild = leftChild
        #keys[parentNode] = GridKey{T,K}(keys[parentNode].index, keys[parentNode].morton_code, keys[parentNode].min, keys[parentNode].max, pleftchild, keys[parentNode].skip)
        keys[parentNode] = GridKey{T,K}(keys[parentNode].min, keys[parentNode].max, pleftchild, keys[parentNode].skip)
        

        if ranger == spec.leaves_count
            #keys[parentNode].skip = K(0)
            zerskip = K(0)
            #keys[parentNode] = GridKey{T,K}(keys[parentNode].min, keys[parentNode].max, keys[parentNode].left, zerskip)

            #TODO temporarily disabled
            keys[parentNode] = GridKey{T,K}(boundingmin, boundingmax, keys[parentNode].left, zerskip)
            
        else
            r = ranger + K(1)
            if delr < delta(r, pos, spec)
                #keys[parentNode].skip = r
                #keys[parentNode] = GridKey{T,K}(keys[parentNode].min, keys[parentNode].max, keys[parentNode].left, r)
                keys[parentNode] = GridKey{T,K}(boundingmin, boundingmax, keys[parentNode].left, r)
            else
                #keys[parentNode].skip = branch_index(r, spec) 
                r = branch_index(r, spec) 
               # keys[parentNode] = GridKey{T,K}(keys[parentNode].min, keys[parentNode].max, keys[parentNode].left, r)
                keys[parentNode] = GridKey{T,K}(boundingmin, boundingmax, keys[parentNode].left, r)
            end
        end
        #KA.@print("happier")
        # keys[parentNode].min = boundingmin
        # keys[parentNode].max = boundingmax

        #keys[parentNode] = GridKey{T,K}(keys[parentNode].index, keys[parentNode].morton_code, boundingmin, boundingmax, keys[parentNode].left, keys[parentNode].skip)
        #keys[parentNode] = GridKey{T,K}(boundingmin, boundingmax, keys[parentNode].left, keys[parentNode].skip)
        i = branch_index(q, spec)


        if i == branch_index(K(1), spec)
            break
        end
        #KA.@print("happiest")
    end

    #min = Adapt.adapt_structureSVector{3, T}(0.0, 0.0, 0.0)
    #max = SVector{3, T}(1.0, 1.0, 1.0)
    #keys[branch_index(1, spec)] = GridKey{T,K}(min, max, keys[branch_index(1, spec)].left, keys[branch_index(1, spec)].skip)

    KA.@synchronize

end

###### Phase 3: traversal


function NaiveDynamics.gpuproximity_test!(neighbors::AbstractVector{Tuple{K, K, T}},  query::IPointPrimitive{T,K}, subject::IPointPrimitive{T, K}, spec::SpheresBVHSpecs{T, K}, squared_radius, query_index, slice_i) where {T, K}
    #a2 = query.position .^ 2
    #for each in eachindex(subjects)
    a = query.index
    b = subject.index

    KA.@print("$a, $b   ")
        if query.index < subject.index
            #dxyz2 = sum(a2 - 2 .* query.position .* subjects[each].position + (subjects[each].position .^ 2))
            dxyz2 = sum( (query.position - subject.position) .^ 2 )
            
            #not demonstrably faster 3/15/2025
            #dxyz2 = sqeuclidean(positions[query_index].position, positions[leafsatoms].position) 

            # only push! new pairs that are close together 
            if dxyz2 < squared_radius

                

                d2 = sqrt( dxyz2 )
                #push!(neighbors, tuple(positions[query_index].index, positions[leafsatoms].index, dxyz, d2))
                #push!(neighbors, tuple(query.index, subject.index, d2))
                #KA.@print("chill")

                #hrm
                #neighbors[query_index*slice_i] = (query.index, subject.index, d2)
                eighbors[query_index+slice_i] = (query.index, subject.index, d2)
            end

        end
    #end

    #return neighbors
end
function NaiveDynamics.gpuoverlap_test(myKey, myPos, spec::SpheresBVHSpecs{T, K}) where {T, K}
    return all(myKey.min .< myPos.position .< myKey.max)

end



@kernel cpu=false function NaiveDynamics.gpuneighbor_traverse!(list::AbstractVector{Tuple{K, K, T}}, @Const(keys::AbstractVector{GridKey{T,K}}), @Const(positions::AbstractVector{IPointPrimitive{T,K}}), @Const(spec::SpheresBVHSpecs{T, K})) where {T, K}
    # pengwing = @groupsize()
    # KA.@print("$pengwing")
    #KA.print(@groupsize())
    #neighbor_vec = parallel_neighbor_buffer(spec)
    
    squared_radius = (spec.neighbor_distance) ^ 2

    # 1 atom per thread
    query_index = KA.@index(Global, Linear)
    #for chunk in 1:threads
         #for query_index in K(chunk):threads:K(length(positions))

         #TODO is this fxn executed on cpu????
            currentKey = branch_index(1, spec)

            while currentKey != 0 # currentKey is the sentinel, end traversal of the given query
                # does query at all overlap with the volume of currentKey
                myKey = keys[currentKey]
                myPos = positions[query_index]
                overlap = gpuoverlap_test(myKey, myPos, spec)
                # overlap = overlap_test(keys, currentKey, query_index, positions, spec)
                #overlap = sum(keys[currentKey].min .< positions[query_index].position .< keys[currentKey].max)
                # overlap = @btime overlap_test($keys, $currentKey, $query_index, $positions, $spec)
                # overlap = @btime newoverlap_test($keys, $currentKey, $query_index, $positions, $spec)
                # return
                # if positions[currentKey].index == 4
                #     KA.@print("$currentKey  ")
                # end

                if overlap

                    #should be replaced with smth like currentKey < maximum index value for any leaf + 1
                    if myKey.left == 0 # currentKey is a leaf node

                        # i like the new method better, but i cannot for the life of me tell if one is better than the other. 
                        #when i repeat the same eval on the same data a million times, i see about a 10% improvement, but the behavior is def broken
                        # also, i like this new method better. Much prettier.
                        #@inbounds newproximity_test!(neighbor_vec[chunk], query_index, currentKey, positions, spec, squared_radius)
                        low = (currentKey-1) * spec.atomsperleaf + 1
                        high = currentKey * spec.atomsperleaf
                        slice = @view positions[low:high]

                        for each in eachindex(slice)

                            if positions[query_index].index < slice[each].index

                                #dxyz2 = sum(a2 - 2 .* query.position .* subjects[each].position + (subjects[each].position .^ 2))
                                dxyz2 = sum( (positions[query_index].position .- slice[each].position) .^ 2 )
                                
                                #not demonstrably faster 3/15/2025
                                #dxyz2 = sqeuclidean(positions[query_index].position, positions[leafsatoms].position) 
                    
                                # only push! new pairs that are close together 
                                if dxyz2 < squared_radius
                                    
                    
                                    d2 = sqrt( dxyz2 )
                                    #push!(neighbors, tuple(positions[query_index].index, positions[leafsatoms].index, dxyz, d2))
                                    #push!(neighbors, tuple(query.index, subject.index, d2))
                                    #KA.@print("chill")
                                    list[query_index+each] = (positions[query_index].index, slice[each].index, d2)
                                end
                    
                            end
                        end


                        # # no idea which!
                        # #slice_i = KA.@index(Group, Linear)
                        # slice_i = KA.@index(Local, Linear)
                        # if slice_i < high+1
                            
                        #     gpuproximity_test!(list, myPos, positions[slice_i],  spec, squared_radius, query_index, slice_i)
                        # end
                        #proximity_test!(neighbor_vec[chunk], query_index, currentKey, positions, spec)
                        currentKey = myKey.skip
                    else #currentKey is a branch node, traverse to the left
                        currentKey = myKey.left
                    end
                else #query is not contained, can cut off traversal on the 'lefts' sequencef
                    currentKey = myKey.skip
                end

            end
        #end
    #end


    @synchronize
    #probably just return the data to the host rather than doing the thread culling  below


    #return 
    #neighbors2 = neighbor_vec[1]
    # for each in 2:1:threads
    #     append!(neighbor_vec[1], neighbor_vec[each] )
    # end

    #not permittedd in GPU prograsmming
    #return neighbor_vec[1]
end

function prune_neighbors!(neighbors)
    counter = 0
    for i in eachindex(neighbors)
        if neighbors[i][1] == 0

        else
            counter+=1
        end
    end
    return [neighbors[each] for each in 1:counter if neighbors[each][1] > 0]
end

###### Phase 4: put it all together


function NaiveDynamics.gpubvh_neighborlist(backend, position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}

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
    


    #data to gpu happens here


    #gpuspec = SpheresBVHSpecs{T, K}(spec.neighbor_distance, spec.atom_count, spec.leaves_count, spec.branches_count, spec.atomsperleaf)
    #println("pre", length(L), " ", sizeof(L))
    
    gpuspec = Adapt.adapt_structure(backend, spec)
    gpustore = Adapt.adapt_structure(backend, store)
    gpupos = Adapt.adapt_structure(backend, pos)
    gpuL = Adapt.adapt_structure(backend, L)
    #println(" $")


    #KA.@print("Happy")
    kernel! = gpubounding_volume_hierarchy!(backend)
    kernel!(gpuL, gpustore, gpuspec, gpupos, ndrange=gpuspec.leaves_count)
    
    #KA.@print("Happy")
    
    
    #println(gpuL)
    #println()

    
    #L = Vector(gpuL)
    # println("post", length(L), " ", sizeof(L))
    # println(typeof(L))
    # for each in eachindex(L)
    #     #println(each, " ", L[each].left, " ", L[each].skip )
    #     println(L[each])
    # end


    # super arbitrary and i imagine results will go to hell if this does not work. def need to have this as a bounds check or something at top of kernel and end
    #neighbor_max = round(K, (spec.atom_count * (spec.neighbor_distance^2) * 1 * 1), RoundUp)
                                                                # number of atom types
                                                                # number of search ranges
    neighbor_max = spec.atom_count * (spec.atom_count - 1) / 2
    #neighbors = [tuple(K(0), K(0), T(0)) for i in 1:neighbor_max]

    neighbors = [tuple(K(0), K(0), T(0)) for i in 1:neighbor_max]
    #neighbors = Vector{Tuple{K, K, T}}(undef, neighbor_max)
    gpuneighbors = Adapt.adapt_structure(backend, neighbors )


    kernel! = gpuneighbor_traverse!(backend)
    kernel!(gpuneighbors, gpuL, gpupos, gpuspec, ndrange=gpuspec.atom_count)

    #TODO this takes a while, can it get better with fwer allocs?
    #println(gpuneighbors)
    neighborslist = prune_neighbors!(Vector(gpuneighbors))
    L = Vector(gpuL)
    #return nothing
    # return data to host here
    
    #println("happy so far")
    return (pairlist=neighborslist, treedata=(tree=L, mortonpositions=pos, quantizedposition=quantized_xyz, store=store))
end






































function NaiveDynamics.gpubvh_neighborlist!(neighborlist, treedata, spec, backend) 

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

    # send data to gpu here

    # our two kernels
    gpubounding_volume_hierarchy!(L, store, spec, pos)

    neighborlist = gpuneighbor_traverse(L, pos, spec)


    return tuple(neighborlist, tuple(L, pos, quantized_xyz, store))
end

end #module