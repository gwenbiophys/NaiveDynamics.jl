export
    SpheresBVHSpecs,
    AxisAlignedBoundingBox,
    AABB,
    QuantizedAABB,
    GridKey,
    update_bvh!,
    build_bvh,
    traverse_bvh

# how to build towards an API that makes it easy to extend a BVH procedure to different kinds of shapes and considering different distances, as in the Noneuclidean paper?
struct SpheresBVHSpecs{T} <: SimulationSpecification
    critical_distance::T
    atoms_count::Int64
    bins_count::Int64

end

function SpheresBVHSpecs(; floattype, interaction_distance, atoms_count, bins_count=atoms_count )
    #TODO incorporate this into the instantiation of the specs type to reduce argument stress throughout
    #if T==Float32
     #   morton_type = Int32
      #  morton_length = [Int32(i) for i in 1:30]
       # morton_string = "000000000000000000000000000000"
        #morton_power = 10
    #elseif T==Float64
     #   morton_type = Int64
      #  mort_length = [Int32(i) for i in 1:63]
       # morton_string = "000000000000000000000000000000000000000000000000000000000000000"
        #morton_power = 21
    #end
    return SpheresBVHSpecs{floattype}(interaction_distance, atoms_count, bins_count)
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

mutable struct GridKey{T} <: AABBGridKey
    index::T
    morton_code::T
end




function generate_aabb(position::Vec3D{T}, spec::SpheresBVHSpecs{T}) where T
    aabb_array = [AABB{T}(i, position[i], position[i] .- spec.critical_distance, position[i] .+ spec.critical_distance) for i in eachindex(position)]


    return aabb_array

end
function update_aabb!(position::Vec3D{T}, spec::SpheresBVHSpecs{T}, aabb_array::Vector{AABB{T}}) where T
    for i in eachindex(aabb_array)
        aabb_array[i].centroid .= position[i]
        aabb_array[i].min .= position[i] .- spec.critical_distance
        aabb_array[i].max .= position[i] .+ spec.critical_distance
    end
    
    return aabb_array
end




# TODO fix this function to only evaluate for the center. forget the boundaries we already know their grid-space values once we know the center's values
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

    return quantized_aabbarray

end

function update_mortoncodes!(L, quantized_aabbs, morton_length, morton_type) 
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
        #println("morton code for atom ", L[each].index, " ", (bitstring(L[each].morton_code)))
    end
end
function assign_mortoncodes(aabb_array, spec::SpheresBVHSpecs{T}, clct) where T
    morton_length = 0::Int64
    morton_string = " "::String
    morton_type = Int32

    #idk how this should be better done to prevent runtime evaluation of a stupid if statement
    if T==Float32

        morton_length = [Int32(i) for i in 1:30]
        morton_string = "000000000000000000000000000000"
        morton_power = 10
    elseif T==Float64
        morton_type = Int64
        mort_length = [Int32(i) for i in 1:63]
        morton_string = "000000000000000000000000000000000000000000000000000000000000000"
        morton_power = 21
    end

    quantized_aabbs = [QuantizedAABB{morton_type}(i, MVector{3, morton_type}(0, 0, 0),  MVector{3, morton_type}(0, 0, 0), MVector{3, morton_type}(0, 0, 0)) for i in 1:spec.bins_count]

    update_gridkeys!(quantized_aabbs, aabb_array, spec, morton_type, clct)

    L = [GridKey{morton_type}(quantized_aabbs[i].index, 0) for i in 1:spec.atoms_count]

    update_mortoncodes!(L, quantized_aabbs, morton_length, morton_type)

    return L

end



function sort_mortoncodes!(L::Vector{GridKey{T}}) where T

    sort!(L, by = x -> x.morton_code)

    return L
end

function create_mortoncodes(position, spec::SpheresBVHSpecs{T}, clct::GenericRandomCollector{T}) where T
    #TODO is the kind of function scoping we want? -- ask again in the refactor
    aabb_array = generate_aabb(position, spec)
    L = assign_mortoncodes(aabb_array, spec, clct)


    sort_mortoncodes!(L)
    
    return L

end

###### Phase 2: Binary radix tree / bvh construction and updating


function commonprefix_length()
end



function build_bvh(position::Vec3D{T}, spec::SpheresBVHSpecs{T}, clct::GenericRandomCollector{T}) where T
    #L is an array of leaf nodes, 1 leaf per atom or 'primitive'
    L = create_mortoncodes(position, spec, clct)
    println(L)
end


function update_bvh!(L, I, position, spec::SpheresBVHSpecs{T}, clct::Collector, aabb_array) where T
    update_aabb!(position, spec, aabb_array)


end


###### Phase 3: Tree traversal

function query_points()

end
function traverse_bvh()

end