module NaiveSIMD

using NaiveDynamics
using StaticArrays
using Polyester
using SIMD


# simd data types are actually larger than SVector, so no thanks don't want it. stick with SVectors and convert when requested
struct SIMDPointPrimitive{T, K}
    index::Vector{K}
    morton_code::Vector{K} #TODO update naming so this becomes one word, maybe even just morton?
    position::Vector{SVector{3, T}}
end


function simdTreeData(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}


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
    
    append!(L, I)


    exptbounding_volume_hierarchy!(L, store, spec, pos)

    #vpos = [SIMDPointPrimitive{T, K}(pos[each].index, pos[each].morton_code, Vec{3,T}((pos[each].position[1], pos[each].position[2], pos[each].position[3]))) for each in eachindex(pos)]
    #vpos = SIMDPointPrimitive{T, K}([pos[each].index for each in eachindex(pos)], [pos[each].morton_code for each in eachindex(pos)], [SVector{3,T}((pos[each].position[1], pos[each].position[2], pos[each].position[3])) for each in eachindex(pos)]) 
    return tuple(L, pos, quantized_xyz, store)
end

Base.@propagate_inbounds function simdproximity_test!(neighbors::Vector{Tuple{K, K, T}}, d2, query, queryi, subjects, subjectsi, spec::SpheresBVHSpecs{T, K}, squared_radius) where {T, K}

    # for each in eachindex(subjects)


    #     if query.index < subjects[each].index

    #         dxyz2 = sum( (query.position - subjects[each].position) ^ 2 )

    #         #not demonstrably faster 3/15/2025
    #         #dxyz2 = sqeuclidean(positions[query_index].position, positions[leafsatoms].position) 

    #         # only push! new pairs that are close together 
    #         if dxyz2 < squared_radius

                

    #             d2 = sqrt( dxyz2 )
    #             #push!(neighbors, tuple(positions[query_index].index, positions[leafsatoms].index, dxyz, d2))
    #             push!(neighbors, tuple(query.index, subjects[each].index, d2))
    #         end

    #     end
    # end
    # d2 = @SVector zeros(spec.atomsperleaf)
    # d2 = zeros(spec.atomsperleaf)
    
    for each in eachindex(subjects)


        #if query.index < subjects[each].index
        if queryi < subjectsi[each]
            #d2 .= sqrt.(sum.( (query .- subjects) .^ 2 ))
            d2[each] = sum( (query - subjects[each]) .^ 2 )


            # only push! new pairs that are close together 
            #if dxyz2 < squared_radius
            if d2[each] < squared_radius
                d2[each] = sqrt(d2[each])
                push!(neighbors, tuple(queryi, subjectsi[each], d2[each]))
            end
                

                #d2 = sqrt( dxyz2 )
                #push!(neighbors, tuple(positions[query_index].index, positions[leafsatoms].index, dxyz, d2))
                #push!(neighbors, tuple(query.index, subjects[each].index, d2))
            #end

        end
    end
    for each in eachindex(d2)
        d2[each] = typemax(T)
    end

    # for each in eachindex(d2)
    #     if queryi < subjectsi[each]
    #         if d2[each] < spec.neighbor_distance
    #             push!(neighbors, tuple(queryi, subjectsi[each], d2[each]))
    #         end
    #     end
    # end




    return neighbors
end



Base.@propagate_inbounds function scalarproximity_test!(neighbors::Vector{Tuple{K, K, T}}, query, queryi, subjects, subjectsi, spec::SpheresBVHSpecs{T, K}, squared_radius) where {T, K}
    #a2 = query.position .^ 2
    for each in eachindex(subjects)


        if queryi < subjectsi[each]
            #dxyz2 = sum(a2 - 2 .* query.position .* subjects[each].position + (subjects[each].position .^ 2))
            dxyz2 = sum( (query - subjects[each]) .^ 2 )

            #not demonstrably faster 3/15/2025
            #dxyz2 = sqeuclidean(positions[query_index].position, positions[leafsatoms].position) 

            # only push! new pairs that are close together 
            if dxyz2 < squared_radius

                

                d2 = sqrt( dxyz2 )
                #push!(neighbors, tuple(positions[query_index].index, positions[leafsatoms].index, dxyz, d2))
                push!(neighbors, tuple(queryi, subjectsi[each], d2))
            end

        end
    end

    return neighbors
end

# as clever as this felt, we spend MUCH more time in assembly and conversion than actually doing work.
Base.@propagate_inbounds function NaiveDynamics.simd_overlaptest(min::SVector{3, T}, pos::SVector{3, T}, maxi::SVector{3, T}) where T

    subject = Vec{8, Float32}((min[1], min[2], min[3], pos[1], pos[2], pos[3],  Float32(0.0), Float32(0.0)))
    query = Vec{8, Float32}((pos[1], pos[2], pos[3], maxi[1], maxi[2], maxi[3],  Float32(1.0), Float32(1.0)))

    return all(subject < query)

end

# this is just much slower than packing in with data to meet 128 bit vec
Base.@propagate_inbounds function simd3_overlaptest(keys, currentKey, query_index, positions, spec::SpheresBVHSpecs{T, K}) where {T, K}
    min = keys[currentKey].min
    pos = positions[query_index].position
    maxi = keys[currentKey].max

    #min = Vec{3, T}((keys[currentKey].min[1], keys[currentKey].min[2], keys[currentKey].min[3]))
    #vpos = Vec{3, T}((positions[query_index].position[1], positions[query_index].position[2], positions[query_index].position[3]))
    #vmaxi = Vec{3, T}((keys[currentKey].max[1], keys[currentKey].max[2], keys[currentKey].max[3]))
    vmin = Vec{3, T}((min[1], min[2], min[3]))
    vpos = Vec{3, T}((pos[1], pos[2], pos[3]))
    vmaxi = Vec{3, T}((maxi[1], maxi[2], maxi[3]))
    return all( (vmin < vpos) & (vpos < vmaxi))

    #return all(subject < query)

end

Base.@propagate_inbounds function simd4_overlaptest(keys, currentKey, query_index, positions, spec::SpheresBVHSpecs{T, K}) where {T, K}


    min = keys[currentKey].min
    pos = positions[query_index].position
    maxi = keys[currentKey].max
    vmin = Vec{4, T}((min[1], min[2], min[3], min[3]))
    vpos = Vec{4, T}((pos[1], pos[2], pos[3], pos[3]))
    vmaxi = Vec{4, T}((maxi[1], maxi[2], maxi[3], maxi[3]))

    return all( (vmin < vpos) & (vpos < vmaxi))

    # min = 
    # pos = 
    # maxi = 
    # vmin = Vec{4, T}((keys[currentKey].min[1], keys[currentKey].min[2], keys[currentKey].min[3], keys[currentKey].min[3]))
    # vpos = Vec{4, T}((positions[query_index].position[1], positions[query_index].position[2], positions[query_index].position[3], positions[query_index].position[3]))
    # vmaxi = Vec{4, T}((keys[currentKey].max[1], keys[currentKey].max[2], keys[currentKey].max[3], keys[currentKey].max[3]))

    # return all( (Vec{4, T}((keys[currentKey].min[1], keys[currentKey].min[2], keys[currentKey].min[3], keys[currentKey].min[3])) 
    #             < Vec{4, T}((positions[query_index].position[1], positions[query_index].position[2], positions[query_index].position[3], positions[query_index].position[3]))) 
    #             & (Vec{4, T}((positions[query_index].position[1], positions[query_index].position[2], positions[query_index].position[3], positions[query_index].position[3])) 
    #             < Vec{4, T}((keys[currentKey].max[1], keys[currentKey].max[2], keys[currentKey].max[3], keys[currentKey].max[3])))
    # ) 

    #return result
    #return all(subject < query)

end

Base.@propagate_inbounds function NaiveDynamics.simd_neighbor_traverse(keys::Vector{GridKey{T,K}}, positions::Vector{IPointPrimitive{T,K}}, spec::SpheresBVHSpecs{T, K}) where {T, K}

    #neighbor_vec = parallel_neighbor_buffer(spec)
    
    threads = K(Threads.nthreads())

    #TODO isn't this structure extremely unfriendly to growing the individual elements at different times to different extents?
    #and what prevents death by thread racing??? certainly, nothing that i have consciously wrote
    # ----- not sure, having this as a vector of references seemed to diminish performance
    neighbor_vec = [Vector{Tuple{K, K, T}}(undef, 0) for i in 1:threads]
    #neighbor_vec = parallel_neighbor_buffer(spec)
    #neighbor_vec = Vector{Vector{Tuple{K, K, T}}}(undef, Threads.nthreads())
    squared_radius = (spec.neighbor_distance) ^ 2


    @batch for chunk in 1:threads
        # d2 = ones(T, spec.atomsperleaf)
        # for each in eachindex(d2)
        #     d2[each] = typemax(T)
        # end
         for query_index in K(chunk):threads:K(length(positions))
            currentKey = branch_index(1, spec)

            while currentKey != 0 # currentKey is the sentinel, end traversal of the given query
                # does query at all overlap with the volume of currentKey
                #overlap = @inbounds simd_overlaptest(keys[currentKey].min, positions[query_index].position, keys[currentKey].max)
                #overlap = @inbounds nextsimd_overlaptest(keys[currentKey].min, positions[query_index].position, keys[currentKey].max)
                # min = keys[currentKey].min
                # pos = positions.position[query_index]
                # maxi = keys[currentKey].max
                # vmin = Vec{4, T}((min[1], min[2], min[3], min[3]))
                # vpos = Vec{4, T}((pos[1], pos[2], pos[3], pos[3]))
                # vmaxi = Vec{4, T}((maxi[1], maxi[2], maxi[3], maxi[3]))

                #overlap = @inbounds presimd_overlaptest(vmin, vpos, vmaxi)
                overlap = @inbounds simd4_overlaptest(keys, currentKey, query_index, positions, spec)
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
                        
                        @inbounds newnewproximity_test!(neighbor_vec[chunk], positions[query_index], slice,  spec, squared_radius)
                        
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

    #return reduce(vcat, neighbor_vec)

    #return neighbors
end

Base.@propagate_inbounds function testsimd_neighbor_traverse(keys::Vector{GridKey{T,K}}, positions::Vector{IPointPrimitive{T,K}}, spec::SpheresBVHSpecs{T, K}) where {T, K}

    #neighbor_vec = parallel_neighbor_buffer(spec)
    
    threads = K(Threads.nthreads())

    #TODO isn't this structure extremely unfriendly to growing the individual elements at different times to different extents?
    #and what prevents death by thread racing??? certainly, nothing that i have consciously wrote
    # ----- not sure, having this as a vector of references seemed to diminish performance
    neighbor_vec = [Vector{Tuple{K, K, T}}(undef, 0) for i in 1:threads]
    #neighbor_vec = parallel_neighbor_buffer(spec)
    #neighbor_vec = Vector{Vector{Tuple{K, K, T}}}(undef, Threads.nthreads())
    squared_radius = (spec.neighbor_distance) ^ 2



    @batch for chunk in 1:threads
        # #d2 = ones(T, spec.atomsperleaf)
        # for each in eachindex(d2)
        #     d2[each] = typemax(T)
        # end
         for query_index in K(chunk):threads:K(length(positions))
            currentKey = branch_index(1, spec)

            while currentKey != 0 # currentKey is the sentinel, end traversal of the given query
                # does query at all overlap with the volume of currentKey
                #overlap = @inbounds simd_overlaptest(keys[currentKey].min, positions[query_index].position, keys[currentKey].max)
                #overlap = @inbounds nextsimd_overlaptest(keys[currentKey].min, positions[query_index].position, keys[currentKey].max)

                overlap = @inbounds simd4_overlaptest(keys, currentKey, query_index, positions, spec)
                #overlap = @inbounds intsimd_overlaptest(keys, currentKey, query_index, positions, spec)
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
                        
                        @inbounds newnewproximity_test!(neighbor_vec[chunk], positions[query_index], slice,  spec, squared_radius)
                        
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

    #return reduce(vcat, neighbor_vec)

    #return neighbors
end

function NaiveDynamics.simdbuild_traverse_bvh(position, spec::SpheresBVHSpecs{T, K}) where {T, K}
    treeData = simdTreeData(position, spec)
    #return expt_neighbor_traverse(treeData.tree, treeData.position, spec)
    return @inbounds simd_neighbor_traverse(treeData[1], treeData[2], spec)
end

function NaiveDynamics.testsimdbuild_traverse_bvh(position, spec::SpheresBVHSpecs{T, K}) where {T, K}
    treeData = simdTreeData(position, spec)
    #return expt_neighbor_traverse(treeData.tree, treeData.position, spec)
    return @inbounds testsimd_neighbor_traverse(treeData[1], treeData[2], spec)
end


end #module
