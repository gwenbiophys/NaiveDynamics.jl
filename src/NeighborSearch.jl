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


struct SpheresBVHSpecs{T, K} <: SimulationSpecification
    critical_distance::T
    atoms_count::Int64
    bins_count::Int64
    morton_length::Vector{K}

end
"""
    function SpheresBVHSpecs(; floattype, interaction_distance, atoms_count, bins_count=atoms_count )

Instantiate a specification towards a BVH of sphere primitives. The ```interaction_distance``` is the maximum interaction distance for the set of 
pairwise interactions that a BVH+traversal algorithm finds the neighbors of. The ```bins_count``` is the number of chunks each axis will be divided by
in order convert particles from real space to grid space. By default, ```bins_count = atoms count```. Morton_length. Though Howard et al., 2019 chose 1023 bins to fit each 
grid axis within a UInt8, instead of here where the integer that fits a grid space axis has the same number of bits as the ```floattype```. 

"""
function SpheresBVHSpecs(; floattype, interaction_distance, atoms_count, bins_count=atoms_count )

    if floattype==Float32
        morton_type = Int32
        morton_length = [Int32(i) for i in 1:30]

        #morton_power = 10
        return SpheresBVHSpecs{floattype, morton_type}(interaction_distance, atoms_count, bins_count, morton_length)
    elseif floattype==Float64
        morton_type = Int64
        mort_length = [Int64(i) for i in 1:63]

        #morton_power = 21
        return SpheresBVHSpecs{floattype, morton_type}(interaction_distance, atoms_count, bins_count, morton_length)
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

# does this need to be mutable?
mutable struct GridKey{T} <: AABBGridKey
    index::T
    morton_code::T
end

# struct internal node
    #leaf indices::Vector(Tuple{T, T})
    #left traversal pointer::Vector{Ref{?}}
    #right traversal pointer::Vector{Ref{?}}

    #leaf indices::Tuple{T, T}
    #left::Ref{Union{internalnode, gridkey}}
    #left::Ref{Union{internalnode, gridkey}}
    #but then this cannot be an array
abstract type NaiveNode end
mutable struct INode{T} <: NaiveNode
    leaf_indices::Tuple{T, T}
    left::Ref{Union{GridKey, NaiveNode}}
    right::Ref{Union{GridKey, NaiveNode}}
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
        #println("morton code for atom ", L[each].index, " ", (bitstring(L[each].morton_code)))
    end
end
function assign_mortoncodes(aabb_array, spec::SpheresBVHSpecs{T, K}, clct) where {T, K}
    morton_length = 0::Int64
    morton_string = " "::String
    #morton_type = K
    #idk how this should be better done to prevent runtime evaluation of a stupid if statement


    quantized_aabbs = [QuantizedAABB{K}(i, MVector{3, K}(0, 0, 0),  MVector{3, K}(0, 0, 0), MVector{3, K}(0, 0, 0)) for i in 1:spec.bins_count]

    update_gridkeys!(quantized_aabbs, aabb_array, spec, K, clct)

    L = [GridKey{K}(quantized_aabbs[i].index, 0) for i in 1:spec.atoms_count]

    update_mortoncodes!(L, quantized_aabbs, spec.morton_length, K)

    return L

end



function sort_mortoncodes!(L::Vector{GridKey{T}}) where T

    sort!(L, by = x -> bitstring(x.morton_code)) 
    #sort!(L, by = x -> string(x.morton_code)) # somtimes does not sort lexicographically

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

"""
    δ(L::Vector{GridKey{K}}, i, j, spec::SpheresBVHSpecs{T,K}) where {T, K}

Compare morton codes in L[i] and L[j] and return the number of common prefix bits, 
wherein the first or second bits of the prefix would be the final or final 2 bits one would increment. 


"""
function δ(i, j, L::Vector{GridKey{K}}, spec::SpheresBVHSpecs{T,K}) where {T, K}# only use where if you are going to use both of these!
    common_prefix = 0
    #println("i in δ ", i)
    #println("j in δ ", j)

    #bandaid fix
    prefixarray = [Int32(i) for i in 1:32]
    reverse!(prefixarray)
    if j > length(L) || j < 1
        return common_prefix -= 1
    else
        for n in prefixarray

            # this is a really heavy handed method.
            a = L[i].morton_code >> (n-1) & 1 != 0
            b = L[j].morton_code >> (n-1) & 1 != 0
            if a == b # if the nth bit of the morton code at positions i and j are the same, then increment common_prefix by 1
                common_prefix += 1
                #println("we incremented: $common_prefix")
            else
                return common_prefix
            end
        end
    end

    return common_prefix
end

function bvh_solver!(L, I, spec)
    # type security in the return of delta is going to be bad
    #δmin = δ(I[1][1], I[1][2], L, spec)
   # println("startingup")
    #println([L[each].morton_code for each in eachindex(L)])
   # for each in eachindex(L)
        #println(bitstring(L[each].morton_code))
    #end

    #println(δmin)
    # this forloop should be parallelizable to 1 processer for each i in L

    for i in 1+1:length(L)-1
        d = 0
        j = 2
        l = 0
        lmax = 1
        zed = 1
        iter = 1
        γ = 0
        
        s = 0
        left = Ref(L, 1)
        right = Ref(L, 1)




       
        #if j < I[1][2] && j > I[1][1] # is there a better method for this bit of control flow? this many if statements should make parallel execution a pain
        #dplus = δ(i, i+1, L, spec)
        #dminus = δ(i, i-1, L, spec)
        d = sign(δ(i, i+1, L, spec) - δ(i, i-1, L, spec))
        δmin = δ(i, i-d, L, spec)
        #println("del min ", δmin)
        #println("d ", d)
        #println("del- ", δ(i, i-1, L, spec))
        #println("del+ ", δ(i, i+d, L, spec))
            #if δ(i, i-1, L, spec) > δmin || δ(i, i-1, L, spec) == δmin
               # d = 1
           # elseif δ(i, i+1, L, spec) > δmin || δ(i, i+1, L, spec) == δmin
               # d = -1
            #end
        #end



        # Here belongs control logic to resolve issues where the morton codes of 2 different atoms are identical.
            # I have no idea how this conflict is resolved.
        if d == 0
            #println("oh shit")
            #println("i ", bitstring(L[i].morton_code))
            #println("i-", bitstring(L[i-1].morton_code))
            #println("i+", bitstring(L[i+1].morton_code))
            return
        end
        #f = l
        while (I[1].leaf_indices[1] <= i + lmax*d <= I[1].leaf_indices[2]) && (δ(i, i + lmax*d, L, spec) > δmin)
            #f = l
            #println(δ(i, i + lmax*d, L, spec))
            #println("entered here ")
            lmax *= 2

        end
        # i am uncertain if this is necesary
        #lmax ÷= 2
        #println("lmax ", lmax)
        #println("post l ", l)

        # here i want to implement the 'binary search where by each iteration in the for loop,
            #iter = lmax / (2 * current step of the iteration)

        while (iter = lmax ÷ (2 * zed)) >= 1
            #iter = lmax ÷ (2 * zed)
            zed +=1
            if δ(i, i + (l+iter)*d, L, spec) > δmin
                l += iter
            end
        end



        j = i + l*d


        δnode = δ(i, j, L, spec)
        zed = 1
        while (iter = l ÷ (2 * zed)) >= 1
            #iter = lmax ÷ (2 * zed)
            zed +=1
            if δ(i, i + (s+iter)*d, L, spec) > δnode
                s += iter
            end
        end




        γ = i + s*d + min(d, 0)

        # this operation, if implemented, generates pointers for traversal that will point to either the next internal node in the hierarchy
            # or to a leaf node
            # i have no idea how to implement this in julia, haha!
            # i also don't know where the gamma or pointer information would go.
            #also, I be reordered upon the generation of its data, and how?
            # Ref(L, γ) etc is the 'poitner' structure
        if min(i, j) == γ
            #left = Ref(L, γ)[] # the leftpointer to carryout traversal
            left = L[γ]
        else
            #left = Ref(I, γ)[]
            left = I[γ]
        end
        if max(i, j) == γ + 1
            #right = Ref(L, γ+1)[]# right = [Lgamma+1]
            right = L[γ+1]
        else 
            #right = Ref(I, γ+1)[]# right = I[gamma+1]
            right = I[γ+1]
        end

        #println("here is i ", i)
        #println(typeof(left))
        #println(typeof(right))

        I[i].leaf_indices = tuple(i, j)#, left, right)
        I[i].left = left
        I[i].right = right



    end



end



function build_bvh(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T}) where {T, K}
    #L is an array of leaf nodes, 1 leaf per atom or 'primitive'
    L = create_mortoncodes(position, spec, clct)
    #println(L)
    #change type of int here to spec.morton_int, probably with Tuple(spec.morton_int[i, length(L)])

    #I = [tuple(i, length(L), Ref(L, i)[], Ref(L, i)[]) for i in 1:length(L)-1] # -1 from L because we want to have nodes = # atoms - 1 in this construction. Howard et al. chose a fixed 1024-1, but eh
    I = [INode{K}(tuple(i, length(L)), Ref(L, i)[], Ref(L, i)[]) for i in 1:length(L)-1]
    #i = 3
    #j = 4
    #println(bitstring(L[i].morton_code))
    #println(bitstring(L[j].morton_code))
    #δ(L, i, j, spec)
    #println(I)
    #println(I)
    bvh_solver!(L, I, spec)

end


function update_bvh!(L, I, position, spec::SpheresBVHSpecs{T}, clct::Collector, aabb_array) where T
    update_aabb!(position, spec, aabb_array)


end


###### Phase 3: Tree traversal

function query_points()

end
function traverse_bvh()

end