module NaiveSIMD

using NaiveDynamics
using StaticArrays
using Polyester
using SIMD


#29 March

struct SIMDPointPrimitive{T, K}
    index::Vector{K}
    morton_code::Vector{K} #TODO update naming so this becomes one word, maybe even just morton?
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
end
function simdTreeData(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}) where {T, K}
    @assert spec.atomsperleaf == 4
    # pos = [IPointPrimitive{T,K}(i, 0, position[i]) for i in 1:spec.atom_count]

    # quantized_xyz = [SVector{3, K}(0, 0, 0) for i in eachindex(position)] 

    # # quantized_positions!(quantized_xyz, pos, spec)

    # # mortoncodes!(pos, quantized_xyz, spec)
    # mortoncodes!(pos, position, spec)

    # sort_mortoncodes!(pos, spec)
    pos = APointPrimitive{T,K}( [i for i in 1:spec.atom_count], 
    [0 for i in 1:spec.atom_count], 
    [position[i] for i in 1:spec.atom_count]
    )
     
    #quantized_xyz = [SVector{3, K}(0, 0, 0) for i in eachindex(position)] 

    # quantized_positions!(quantized_xyz, pos, spec)

    mortoncodes!(pos, position, spec)

    #sort_mortoncodes!(pos, spec)
    mortsort = sortperm(pos.morton_code)
    permute!(pos.morton_code, mortsort)
    permute!(pos.index, mortsort)
    permute!(pos.position, mortsort)
    #1 alloc per item in this generator expression
    #store = [Threads.Atomic{K}(0) for i in 1:spec.branches_count]
    store = [K(0) for i in 1:spec.branches_count]
    # @time for each in eachindex(store)
    #     store[each] = Threads.Atomic{K}(store[each])
    # end


    L = leafcluster_primitives(pos, spec)




    #I = [GridKey{T, K}(SVector{3, T}(0.0, 0.0, 0.0), SVector{3, K}(0.0, 0.0, 0.0), 0, 0) for i in 1:spec.branches_count]
    
    #append!(L, I)


    bounding_volume_hierarchy!(L, store, spec, pos)

    simdpos = SIMDPointPrimitive{T,K}(  [pos.index[i] for i in 1:spec.atom_count], 
                                        [pos.morton_code[i] for i in 1:spec.atom_count], 
                                        [pos.position[i][1] for i in 1:spec.atom_count],
                                        [pos.position[i][2] for i in 1:spec.atom_count],
                                        [pos.position[i][3] for i in 1:spec.atom_count]
    )
    #return tuple(L, pos, quantized_xyz, store)
    return tuple(L, simdpos, mortsort, store)
end


Base.@propagate_inbounds function gather_A_direct(a)

    x = a[1]
    y = a[2]
    z = a[3]

    return (    Vec{16, Float32}((x[1],x[1],x[1],x[1], x[2],x[2],x[2],x[2], x[3],x[3],x[3],x[3], x[4],x[4],x[4],x[4])),
                Vec{16, Float32}((y[1],y[1],y[1],y[1], y[2],y[2],y[2],y[2], y[3],y[3],y[3],y[3], y[4],y[4],y[4],y[4])),
                Vec{16, Float32}((z[1],z[1],z[1],z[1], z[2],z[2],z[2],z[2], z[3],z[3],z[3],z[3], z[4],z[4],z[4],z[4]))
    )   
end

Base.@propagate_inbounds function gather_B_direct(a)

    x = a[1]
    y = a[2]
    z = a[3]
    return (Vec{16, Float32}((x[1],x[2],x[3],x[4], x[1],x[2],x[3],x[4], x[1],x[2],x[3],x[4], x[1],x[2],x[3],x[4])),
            Vec{16, Float32}((y[1],y[2],y[3],y[4], y[1],y[2],y[3],y[4], y[1],y[2],y[3],y[4], y[1],y[2],y[3],y[4])),
            Vec{16, Float32}((z[1],z[2],z[3],z[4], z[1],z[2],z[3],z[4], z[1],z[2],z[3],z[4], z[1],z[2],z[3],z[4]))
    )

end

Base.@propagate_inbounds function simdonecluster_proximitytest!(list::Vector{Tuple{K, K, T}}, cluster,  spec::SpheresBVHSpecs{T, K}, squared_radius) where {T, K}

    prealloc = Vec{6, Float32}((0, 0, 0, 0, 0, 0))
    a = (1, 1, 1, 2, 2, 3)
    b = (2, 3, 4, 3, 4, 4)
    #Vec6 may be slower than straight scalar
    vec_a = (    Vec{6, T}((cluster[2][1][1], cluster[2][1][1], cluster[2][1][1], cluster[2][1][2], cluster[2][1][2], cluster[2][1][3])),
                Vec{6, T}((cluster[2][2][1], cluster[2][2][1], cluster[2][2][1], cluster[2][2][2], cluster[2][2][2], cluster[2][2][3])),
                Vec{6, T}((cluster[2][3][1], cluster[2][3][1], cluster[2][3][1], cluster[2][3][2], cluster[2][3][2], cluster[2][3][3]))
    )
    vec_b = (    Vec{6, T}((cluster[2][1][2], cluster[2][1][3], cluster[2][1][4], cluster[2][1][3], cluster[2][1][4], cluster[2][1][4])),
                Vec{6, T}((cluster[2][2][2], cluster[2][2][3], cluster[2][2][4], cluster[2][2][3], cluster[2][2][4], cluster[2][2][4])),
                Vec{6, T}((cluster[2][3][2], cluster[2][3][3], cluster[2][3][4], cluster[2][3][3], cluster[2][3][4], cluster[2][3][4]))
    )
    
    #posa = 
    for i in eachindex(vec_a)
        prealloc += (vec_a[i] - vec_b[i]) ^ 2
    end
    better = @inbounds NTuple{6, T}(sqrt(prealloc))
    for i in eachindex(better)
        #println(clusterA[1][indexA[i]], " ",  clusterB[1][indexB[i]], " ",  better[i])
        if better[i] < spec.neighbor_distance
            push!(list, (cluster[1][a[i]], cluster[1][b[i]], better[i]))
        end
    end




    # for i in eachindex(apos)
    #     #dxyz2=T(0.0)
    #     for j in eachindex(a[2])
    #         dxyz2 += (cluster[2][j][a[i]] - cluster[2][j][b[i]]) ^2
    #     end
    #     if dxyz2 < squared_radius
    #         d2 = sqrt(dxyz2)
    #         push!(neighbors, (cluster[1][a[i]], cluster[1][b[i]], d2))
    #     end
    # end

    # for i in 1:spec.atomsperleaf-1
    #     for j in i+1:spec.atomsperleaf
    #         dxyz2 = sum( (cluster[2][i] - cluster[2][j]) .^2 )
    #         if dxyz2 < squared_radius
    #             d2 = sqrt(dxyz2)

    #             # maybe instead of push this makes a cluster of pairings or some comprehension, and then we make and append a vector below??
    #             # intuition says compiler would give better performance IF push asnd maybe ifdxyz2 were not here
    #             push!(neighbors, (cluster[1][i], cluster[1][j], d2))
    #         end
    #     end

    # end
    #append!(neighbors, [tuple(i,j,DistFxn(iterators)) for iter in iterator if dist < threshold])
    return list

end

Base.@propagate_inbounds function simdtwocluster_proximitytest!(list, clusterA, clusterB, prealloc, spec::SpheresBVHSpecs{T, K}) where{T, K}

    vec_a = @inbounds gather_A_direct(clusterA[2])
    vec_b = @inbounds gather_B_direct(clusterB[2])

    preallloc = Vec{16, T}((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    for i in eachindex(vec_a)
        preallloc += (vec_a[i] - vec_b[i]) ^ 2
    end
    #preallloc = sqrt(preallloc)
    #TODO last alloc can just be a locally defined NTuple. maybe NTuple is our direction to append! working well
    better = @inbounds NTuple{16, T}(sqrt(preallloc))
   # @inbounds vstore(prealloc, better, 1)
    counter = 0

#    there should just be another fxn called gather indices, and the prior can be gather_simd_positions and gather_positions
    indexA = (1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4,4)
    indexB = (1,2,3,4, 1,2,3,4, 1,2,3,4, 1,2,3,4)
    
    for i in eachindex(better)
        #println(clusterA[1][indexA[i]], " ",  clusterB[1][indexB[i]], " ",  better[i])
        if better[i] < spec.neighbor_distance
            push!(list, (clusterA[1][indexA[i]], clusterB[1][indexB[i]], better[i]))
        end
    end
    # for i in eachindex(clusterA[1])
    #     for j in eachindex(clusterB[1])
    #         counter+=1
    #         if lastalloc[counter] < spec.neighbor_distance #threshold
    #             push!(list, (clusterA[1][i], clusterB[1][j], lastalloc[counter]))
    #         end
    #     end
    # end

    #prealloc = Vec{16, T}((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

    # it just,, kind of works
    #lastalloc .= 0.0
    return list
end

Base.@propagate_inbounds function simdneighbor_traverse(keys::Vector{GridKey{T,K}}, positions::SIMDPointPrimitive{T,K}, spec::SpheresBVHSpecs{T, K}) where {T, K}

    #neighbor_vec = parallel_neighbor_buffer(spec)
    
    threads = K(Threads.nthreads())

    #TODO isn't this structure extremely unfriendly to growing the individual elements at different times to different extents?
    #and what prevents death by thread racing??? certainly, nothing that i have consciously wrote
    # ----- not sure, having this as a vector of references seemed to diminish performance
    #neighbor_vec = [Vector{Tuple{K, K, SVector{3, T}, T}}(undef, 0) for i in 1:threads]
    neighbor_vec = [Vector{Tuple{K, K, T}}(undef, 0) for i in 1:threads]#parallel_neighbor_buffer(spec)
    #neighbor_vec = Vector{Vector{Tuple{K, K, T}}}(undef, Threads.nthreads())
    squared_radius = (spec.neighbor_distance) ^ 2
    prior_leaf_start = 0
    # b = @views (positions.index[1:5], positions.position[1:5])
    # println(typeof(b))

#TODO make query_ into inquisitor, and target_ into quarry_ ??? lol
#difficult because query_node is inquired upon, and becomes and inquisitor unto the following nodes
    @batch for chunk in 1:threads
        #thread local intermediates for twocluster
        vec_store = Vec{16, T}((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
        #result_store = zeros(T, 16)
         for query_index in K(chunk):threads:K(spec.leaves_count)
            ## the query is a fixed identity that will always be a leaf
            query_leaf = keys[query_index]
            
            low = (query_index-1) * spec.atomsperleaf + 1
            high = query_index * spec.atomsperleaf
            #persistent set of pointprimitives for duration of traversal with query_leaf
            #what if this were a staticarray?? huh, huh huhhhhhhh. dumb ideas
            query_cluster = @views ( positions.index[low:high],  (positions.x[low:high], positions.y[low:high], positions.z[low:high]))
            
            ## the target is a changing identity that is either sentinel, internal, or leaf node
            target_index = query_leaf.skip


            @inbounds simdonecluster_proximitytest!(neighbor_vec[chunk], query_cluster, spec, squared_radius)
            # well, this did not quuuite work
            # query_cluster_a = @view positions[low:high-1]
            # query_cluster_b = @view positions[low+1:high]
            # @inbounds twocluster_proximitytest!(neighbor_vec[chunk], query_cluster_a, query_cluster_b,  spec, squared_radius)

            while target_index != 0 # currentKey is the sentinel, end traversal of the given query
                # does query at all overlap with the volume of currentKey
                # myKey = keys[currentKey]
                # myPos = positions[query_index]
                target_node = keys[target_index]

                overlap = @inbounds aabb_overlap_test(query_leaf, target_node, spec)
                # overlap = overlap_test(keys, currentKey, query_index, positions, spec)
                #overlap = sum(keys[currentKey].min .< positions[query_index].position .< keys[currentKey].max)
                # overlap = @btime overlap_test($keys, $currentKey, $query_index, $positions, $spec)
                # overlap = @btime newoverlap_test($keys, $currentKey, $query_index, $positions, $spec)
                # return

                if overlap

                    if target_node.left == 0 # currentKey is a leaf node

                        low = (target_index-1) * spec.atomsperleaf + 1
                        high = target_index * spec.atomsperleaf
                        target_cluster = @views ( positions.index[low:high],   (positions.x[low:high], positions.y[low:high], positions.z[low:high]))
                        
                        @inbounds simdtwocluster_proximitytest!(neighbor_vec[chunk], query_cluster, target_cluster, vec_store, spec)
                        #return
                        #proximity_test!(neighbor_vec[chunk], query_index, currentKey, positions, spec)
                        target_index = target_node.skip
                    else #currentKey is a branch node, traverse to the left
                        target_index = target_node.left
                    end
                else #query is not contained, can cut off traversal on this 'left' section of the tree
                    target_index = target_node.skip
                end

            end
        end
    end
    return reduce(vcat, neighbor_vec)
end

function NaiveDynamics.simdbuild_traverse_bvh(position, spec::SpheresBVHSpecs{T, K}) where {T, K}
    treeData = simdTreeData(position, spec)

    return @inbounds simdneighbor_traverse(treeData[1], treeData[2], spec)
end

end #module
