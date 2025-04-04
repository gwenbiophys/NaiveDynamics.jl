#export
    # SpheresBVHSpecs,
    # AxisAlignedBoundingBox,
    # AABB,
    # QuantizedAABB,
    # GridKey,
    # update_bvh!,
    # build_bvh,
    # traverse_bvh
    # #Atomic

# how to build towards an API that makes it easy to extend a BVH procedure to different kinds 
# of shapes and considering different distances, as in the Noneuclidean paper?


# struct SpheresBVHSpecs{T, K} <: SimulationSpecification
#     critical_distance::T
#     atoms_count::Int64
#     bins_count::Int64
#     morton_length::K

# end

function SpheresBVHSpecs(; floattype, interaction_distance, atoms_count, bins_count=atoms_count )

    if floattype==Float32
        morton_type = Int32
        morton_length = Int32(10)

        #morton_power = 10
        return SpheresBVHSpecs{floattype, morton_type}(interaction_distance, atoms_count, bins_count, morton_length)
    elseif floattype==Float64
        morton_type = Int64
        mort_length = Int64(21)

        #morton_power = 21
        return SpheresBVHSpecs{floattype, morton_type}(interaction_distance, atoms_count, bins_count, morton_length)
    end
    
end

##### Phase 1: morton code construction and updating


# abstract type AxisAlignedBoundingBox end 
# struct AABB{T} <: AxisAlignedBoundingBox where T
#     index::Int64
#     centroid::MVector{3, T}
#     min::MVector{3, T}
#     max::MVector{3, T}
# end

# abstract type AABBGridKey end
# struct QuantizedAABB{T} <: AxisAlignedBoundingBox
#     index::T
#     #this would be better, but it's easier to make something work with MVecs
#     #maybe if the rest of this magically turns out to be realllly easy, i'll work here some more
#     #min::Tuple{T, T, T}
#     #max::Tuple{T, T, T}
#     centroid::MVector{3, T}

#     min::MVector{3, T}
#     max::MVector{3, T}
# end

# # does this need to be mutable?
# mutable struct GridKey{T, K} <: AABBGridKey
#     index::T
#     morton_code::T
#     min::MVector{3, K}
#     max::MVector{3, K}
#     parent_INode::T # this may be repurposed as the index of the skip connection
#     #left::Ref{Union{GridKey, NaiveNode}}
#     #right::Ref{Union{GridKey, NaiveNode}}
#     left::T
#     skip::T
# end

# # struct internal node
#     #leaf indices::Vector(Tuple{T, T})
#     #left traversal pointer::Vector{Ref{?}}
#     #right traversal pointer::Vector{Ref{?}}

#     #leaf indices::Tuple{T, T}
#     #left::Ref{Union{internalnode, gridkey}}
#     #left::Ref{Union{internalnode, gridkey}}
#     #but then this cannot be an array
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

#mutable struct Atomic{T}; @atomic x::T; end




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

function update_mortoncodes2!(L, quantized_aabbs, morton_length, morton_type) 

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
            for dim in t3:-1:t1 #iterate backwards

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


    quantized_aabbs = [QuantizedAABB{K}(i, MVector{3, K}(0, 0, 0),  MVector{3, K}(0, 0, 0), MVector{3, K}(0, 0, 0)) for i in 1:spec.bins_count]::Vector{QuantizedAABB{K}}

    update_gridkeys!(quantized_aabbs, aabb_array, spec, K, clct)

    L = [GridKey{K, T}(quantized_aabbs[i].index, 0, aabb_array[i].min, aabb_array[i].max, 0, 0, 0) for i in 1:spec.atoms_count]

    update_mortoncodes2!(L, quantized_aabbs, spec.morton_length, K)


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

    sort!(L, by = x -> bitstring(x.morton_code)) # sorts lexicographically both the binary and the integer
    #sort!(L, by=x -> count(c -> c == '1', bitstring(x.morton_code)))
    #sort!(L, by = x -> x.morton_code), alg=RadixSort #wont run, RadixSort does not have iterate defined

    return L
end

function create_mortoncodes(position, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T}) where {T, K}
    #TODO is the kind of function scoping we want? -- ask again in the refactor
    aabb_array = generate_aabb(position, spec)
    L = assign_mortoncodes(aabb_array, spec, clct)



    sort_mortoncodes!(L)




    
    return L

end

###### Phase 3: Binary radix tree based on Karras 2012 and Howard 2016 and 2019

"""
    δ(L::Vector{GridKey{K}}, i, j, spec::SpheresBVHSpecs{T,K}) where {T, K}

Compare morton codes in L[i] and L[j] and return the number of common prefix bits, 
wherein the first or second bits of the prefix would be the final or final 2 bits one would increment. 
This is the Karras 2012 method.


"""
function δ(i, j, L::Vector{GridKey{K, T}}, spec::SpheresBVHSpecs{T,K}) where {T, K}# only use where if you are going to use both of these!
    common_prefix = 0
    #println("i in δ ", i)
    #println("j in δ ", j)

    #bandaid fix
    prefixarray = [Int32(i) for i in 1:32]
    reverse!(prefixarray)
    if j > length(L) || j < 1
        return common_prefix = -1
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





function solve_node_1(L, I, spec)

    d = 1
    i = 1
    j = length(L)
    l = 0
    lmax = 1
    zed = 1
    iter = 1
    γ = 0
    
    s = 0
    left = Ref(L, 1)
    right = Ref(L, 1)



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
        L[γ].parent_INode = i 
        #I[i].visits += 1
    else
        #left = Ref(I, γ)[]
        left = I[γ]
        I[γ].parent_INode = i
    end
    if max(i, j) == γ + 1
        #right = Ref(L, γ+1)[]# right = [Lgamma+1]
        right = L[γ+1]
        L[γ+1].parent_INode = i # i am uncertain if it correctly follows convention to have gamma plus 1 here

    else 
        #right = Ref(I, γ+1)[]# right = I[gamma+1]
        right = I[γ+1]
        I[γ+1].parent_INode = i
    end

    #println("here is i ", i)
    #println(typeof(left))
    #println(typeof(right))

    I[i].leaf_indices = tuple(i, j)#, left, right)
    I[i].left = left
    I[i].right = right
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
    lengthL = length(L)
    #solve_node_1(L, I, spec)

    for i in 1:lengthL-1 # Threads.@threads when the atomic Inode.visits is figured out.
        #println(i)
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
        while (1 <= i + lmax*d <= lengthL) && (δ(i, i + lmax*d, L, spec) > δmin)
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
            L[γ].parent_INode = i 
            #I[i].visits += 1
        else
            #left = Ref(I, γ)[]
            left = I[γ]
            I[γ].parent_INode = i
        end
        if max(i, j) == γ + 1
            #right = Ref(L, γ+1)[]# right = [Lgamma+1]
            right = L[γ+1]
            L[γ+1].parent_INode = i # i am uncertain if it correctly follows convention to have gamma plus 1 here

        else 
            #right = Ref(I, γ+1)[]# right = I[gamma+1]
            right = I[γ+1]
            I[γ+1].parent_INode = i
        end
        #println("here is i ", i)
        #println(typeof(left))
        #println(typeof(right))

        I[i].leaf_indices = tuple(i, j)#, left, right)
        I[i].left = left
        I[i].right = right
    end
end
function internal_boundaries(iter, L, I, spec)
    n = iter
    #if 
    p = L[iter].parent_INode # p is the index of the INode that is L[i] 's parent
    #v::Int32 = 0 #maybe store as an array? one for each Inode? represented as I[p].visits
    if p == 0 return end
    # Question TODO what was my goal here?
    if I[p].left == L[n]
       local s = I[p].right
    else
       local s = I[p].left
    end
    #s = # at this moment, i don't know if the nth leaf node is the left or right sibling to the parent

    ######attempt to recreate Karras 2012, Howard 2016
    #this should be an atomic for visits part
    if I[p].visits != 0 #this will never work, need to be reformed
        #merge n and s, whatever the hell that means
        n = p
    end

    ######my attempt
    #need repeated evaluation at each leaf
    while p > 0 # continue climbing up hierarchy until reaching the root, except this prevents the root AABB'ing successfully
        #I[p].visits += 1
        
        left = sum(I[p].left[].min)
        right = sum(I[p].right[].min)
        if 0 < left  && left < right #|| right == 0 # can this be neatly squashed into a ternary? i think so, but with a surrounding if statement or 2
            # if the right value is zero, then we do not! want to set the min value because wwe dont have the sibling comparison, so the else should
            # do nothing, and leave it up to a different 'thread' so bvh_solver
            # we should prune it up by evaluating only the parents of 2 leaves, wiping out the leaf and branch combos
            I[p].min = I[p].left[].min
        elseif 0 < right < left
            I[p].min = I[p].right[].min
        else
            #println("oops in the min node boundaries")
            #println(left)
            #println(right)
        end

        left = sum(I[p].left[].max)
        right = sum(I[p].right[].max)
        if 0 < left < right
            I[p].max = I[p].left[].max
        elseif 0 < right < left
            I[p].max = I[p].right[].max
        else
           # println("oops in the max node boundaries")
           # println(left)
           # println(right)
        end

                        #a = (I[p].left.min < I[p].right.min) * I[p].left.min  
                        #b = (I[p].left.min > I[p].right.min) * I[p].right.min#oh my god, we have to figure out if it is a leaf or inode, which sibling has the min coord and which has max coord. this is so broken
            # i cannot fathom trying to do this this way
        #I[p].max = (I[p].left.max < I[p].right.max) * I[p].left.max + (I[p].left.max > I[p].right.max) * I[p].right.max
        #okay, now how do we get a new uhhhhhhhhhh parent to look at? fuckkkk. 
        #Howard 2016 implies that something that leads us to the new parent is created from mergine the currently considered thread with the sibling
        p = I[p].parent_INode
    end
end

function boundaries_wrapper(L, I, spec)

    #TODO fix these dirty initializations
        # whole thing has to be healed up to allow for initialization and then just update after
        #p::Int32 = 0
        #n::Int32 = 0
        for i in eachindex(L)
            internal_boundaries(i, L, I, spec)
        end

end

function update_bvh!(L, I, position, spec::SpheresBVHSpecs{T}, clct::Collector, aabb_array) where T
    update_aabb!(position, spec, aabb_array)
end

###### Phase 5: Tree traversal

function query_points()

end
function stackless_traverse(L, I, spec::SpheresBVHSpecs{T, K}) where {T, K}
    remaining_nodes = length(L) # meant to be an atomic, proceed until this value is zero.
    neighbors::Vector{Tuple{Int64, Int64, T, T, T, T}} = []

    ### from Howard 2016
    for atom in L # also a parallelizable loop
        n = I[1]
        not_leaf = 1
        while not_leaf > 0
            not_leaf -= 1
            if prod(atom.min .> n.min) || prod(atom.max .< n.max)
                if typeof(n.left[]) != GridKey
                    n = n.left[]
                    not_leaf += 1
                else
                    dx = atom.centroid[1] - n.centroid[1]
                    dy = atom.centroid[2] - n.centroid[2]
                    dz = atom.centroid[3] - n.centroid[3]
                    d2 = sqrt(dx^2 + dy^2 + dz^2)  


                    if d2 <= spec.critical_distance
                        push!(neighbors, tuple(atom.index, n.index, dx, dy, dz, d2))
                    end
                end
            else
                n = n.rope
                not_leaf += 1
            end
       
                # n = I[1] #root node of bvh

        #while remaining_nodes > 0
            # if aabb of atom overlaps aabb of n

            
                #if typeof(n) != GridKey
                   # n = n.left
                #else
                    #for each in 


        end
    end



end

function inner_traverse(neighbors, atom, n, position, L, I, spec::SpheresBVHSpecs{T, K}) where {T, K}

    while typeof(n) != GridKey{K, T} #!== GridKey{T, K}
        #println(typeof(n.left[]))
        #println(typeof(n.right[]))
        if prod(atom.min .> n.left[].min) || prod(atom.max .< n.left[].max)
            #println("here 3")
            if typeof(n.left[]) == GridKey{K, T}
                #println("here after 3")
                dx = position[atom.index][1] - position[n.left[].index][1]
                dy = position[atom.index][2] - position[n.left[].index][2]
                dz = position[atom.index][3] - position[n.left[].index][3]
                d2 = sqrt(dx^2 + dy^2 + dz^2)  
                #println(d2)


                if d2 <= spec.critical_distance
                    #println("here way after 3")
                    push!(neighbors, tuple(atom.index, n.left[].index, dx, dy, dz, d2))
                    #println(neighbors)
                end
                return
            else
                n = n.left[]
            end
        else
            #println("here too")
            if typeof(n.right[]) == GridKey{K, T}
                #println("here after too")
                dx = position[atom.index][1] - position[n.right[].index][1]
                dy = position[atom.index][2] - position[n.right[].index][2]
                dz = position[atom.index][3] - position[n.right[].index][3]
                d2 = sqrt(dx^2 + dy^2 + dz^2)  


                if d2 <= spec.critical_distance
                    push!(neighbors, tuple(atom.index, n.right[].index, dx, dy, dz, d2))
                end
                return
            else
                #println("here")
                n = n.right[]
            end
        end
    end
end

function traverse_bvh1(position, L, I, spec::SpheresBVHSpecs{T, K}) where {T, K}
    remaining_nodes = length(L) # meant to be an atomic, proceed until this value is zero.
    neighbors::Vector{Tuple{Int64, Int64, T, T, T, T}} = []


    for atom in L
        n = I[1]
        ##println(typeof(n), typeof(n) === INode{K, T})
       
        inner_traverse(neighbors, atom, n, position, L, I, spec)

        # for once wwe have ran out of internal nodes, process the last leaf(?)
        if typeof(n) == GridKey{K, T}
            dx = position[atom.index][1] - position[n.index][1]
            dy = position[atom.index][2] - position[n.index][2]
            dz = position[atom.index][3] - position[n.index][3]
            d2 = sqrt(dx^2 + dy^2 + dz^2)  


            if d2 <= spec.critical_distance
                push!(neighbors, tuple(atom.index, n.index, dx, dy, dz, d2))
            end
        end
    end
    return neighbors
end

function traverse_bvh2(L, I, spec::SpheresBVHSpecs{T, K}) where {T, K}
    remaining_nodes = length(L) # meant to be an atomic, proceed until this value is zero.
    neighbors::Vector{Tuple{Int64, Int64, T, T, T, T}} = []


    #### take 2
    for atom in L
        n = I[1]
        while typeof(n) != GridKey
            if prod(atom.min .> n.min) || prod(atom.max .< n.max) # we use 'or' so we can evaluate all cases of any overlap between atom AABB and node AABB
                if typeof(n) != GridKey
                    n = n.left[] # this should not work because we have not tested if it fits more in the left or the right
                else
                    dx = atom.centroid[1] - n.centroid[1]
                    dy = atom.centroid[2] - n.centroid[2]
                    dz = atom.centroid[3] - n.centroid[3]
                    d2 = sqrt(dx^2 + dy^2 + dz^2)  


                    if d2 ⋜ spec.critical_distance
                        push!(neighbors, tuple(atom.index, n.index, dx, dy, dz, d2))
                    end
                end

            else
                if typeof(n.right[]) != GridKey
                    n = n.right[]
                else
                    dx = atom.centroid[1] - n.centroid[1]
                    dy = atom.centroid[2] - n.centroid[2]
                    dz = atom.centroid[3] - n.centroid[3]
                    d2 = sqrt(dx^2 + dy^2 + dz^2)  


                    if d2 ⋜ spec.critical_distance
                        push!(neighbors, tuple(atom.index, n.index, dx, dy, dz, d2))
                    end
                end
            end
        end
    end
end



###### Phase 6: put it all together

# change this name?
function build_bvh(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T}) where {T, K}

    L = create_mortoncodes(position, spec, clct)::Vector{GridKey{K, T}} 

    I= [INode{K, T}(tuple(i, length(L)), 0, 0, 0,  MVector{3, T}(0, 0, 0), MVector{3, T}(0, 0, 0), 0) for i in 1:length(L)-1]::Vector{INode{K, T}} 

end









####### sort perm implementation -- it is legally better than struct sort, but I don't like it as much


function update_quantized_positions(index, quantized, perm, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T}) where {T, K}
    for each in eachindex(quantized.x)
        quantized.x[index.x[perm.x[each]]] = K(each) #TODO i hope this is right, i have no idea
    end
    for each in eachindex(quantized.y)
        quantized.y[index.y[perm.y[each]]] = K(each)
    end
    for each in eachindex(quantized.z)
        quantized.z[index.z[perm.z[each]]] = K(each)
    end

end

function create_mortoncodes(position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T}) where {T, K}
    position_xyz = XYZVectors{T}([position[i][1] for i in 1:spec.leaves_count],
                                 [position[i][2] for i in 1:spec.leaves_count],
                                 [position[i][3] for i in 1:spec.leaves_count]
    )
    index_xyz = XYZVectors{K}([i for i in 1:spec.leaves_count],
                              [i for i in 1:spec.leaves_count],
                              [i for i in 1:spec.leaves_count]
    )
    perm_xyz = XYZVectors{K}(
        sortperm(position_xyz.x),
        sortperm(position_xyz.y),
        sortperm(position_xyz.z)
    )
    #partialsort
    quantized_xyz = XYZVectors{K}(zeros(K, spec.leaves_count),
                              zeros(K, spec.leaves_count),
                              zeros(K, spec.leaves_count)
    )
    update_quantized_positions(index_xyz, quantized_xyz, perm_xyz, spec, clct)

    #TODO this has to be fixed with the permuters and syntax
    L = [GridKey{T, K}(i, 0, 
            #MVector{3, T}(position_xyz.x[perm_xyz.x[i]] - spec.critical_distance, position_xyz.y[perm_xyz.y[i]] - spec.critical_distance, position_xyz.z[perm_xyz.z[i]] - spec.critical_distance),
            #MVector{3, T}(position_xyz.x[perm_xyz.x[i]] + spec.critical_distance, position_xyz.y[perm_xyz.y[i]] + spec.critical_distance, position_xyz.z[perm_xyz.z[i]] + spec.critical_distance),
            SVector{3, T}(position[i] .- spec.critical_distance), SVector{3, T}(position[i] .+ spec.critical_distance),
            0, 0) for i in 1:spec.leaves_count
    ]

    update_mortoncodes!(L, quantized_xyz, spec.morton_length, K)
    #TODO this was greyed out ->sort_mortoncodes!(L, spec)#@time "mortons" sort_mortoncodes!(L)
    store = [Base.Threads.Atomic{Int64}(0) for i in 1:spec.branches_count]
    
    return tuple(Ref(L), Ref(position_xyz), Ref(index_xyz), Ref(perm_xyz), Ref(quantized_xyz), Ref(store))
end




function bcreate_mortoncodes(; position, thresh=0.01, usebtime=false)
    clct = GenericRandomCollector(; floattype=Float32,
                                        objectnumber=length(position),
                                        minDim=tuple(0.0, 0.0, 0.0),
                                        maxDim=tuple(1.0, 1.0, 1.0),
                                        temperature=0.01,
                                        randomvelocity=false,
                                        minmass=1.0,
                                        maxmass=5.0,
                                        minimumdistance=0.0001,
                                        mincharge=-1f-9,
                                        maxcharge=1f-9,
                                        pregeneratedposition=true
    )
    clxn = collect_objects(clct, position=position)
    #TODO greyed out -># bvhspec = SpheresBVHSpecs(; floattype=Float32, 
    #                             critical_distance=thresh, 
    #                             leaves_count=length(position) 
    # )
    atoms = length(position)
    println(atoms, " atoms")
    println("permsort")
    if usebtime
        println("permsort")
        @btime create_mortoncodes($clxn.position, $bvhspec, $clct)
        println("direct sort")
        @btime TreeData($clxn.position, $bvhspec, $clct)
    else
        @time create_mortoncodes(clxn.position, bvhspec, clct)
        println("direct sort")
        @time TreeData(clxn.position, bvhspec, clct)
    end

    @profview for i in 1:4000
        create_mortoncodes(clxn.position, bvhspec, clct)
    end
    @profview for i in 1:4000
        TreeData(clxn.position, bvhspec, clct)
    end

    return nothing
end