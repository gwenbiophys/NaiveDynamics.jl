##### data and tester funcs
# TODO get this to work instead of the recursive traversal pit
function self_traverse(keys, neighbors, position, spec::SpheresBVHSpecs{T, K}) where {T, K}
    currentKey = 1 + spec.leaves_count
    selfVisitedLeaves = [false for each in eachindex(position)]

    while true



        if currentKey == 0 # is sentinel node
            break
        end
    end

    if sum(selfVisitedLeaves) == spec.leaves_count
        return true
    else
        return false
    end


end


#### errata
# those times when you have to traverse down every path
# total sum of stack overflows / runaway Julia before I fixed this function and all its variants
# and caught recursion generating data conditions: 1
# to adequately prevent recursion, we would need a companion array for both the leaves and each of their connections
# as in, evaluate 'has this rope connection been used to update the traversal path already'
# and that's maybe too much work, so let's just run the above counter anyway.
function recursive_traversal(index, keys::Vector{GridKey{T, K}}, wasVisited::Vector{Bool}, spec::SpheresBVHSpecs{T, K}) where {T, K}
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
            if (skip == 0) && (left == 0) # is a sentinel, prepare to exit
                index == skip

            elseif (skip == 0) && (left != 0)# is a right most node, proceed through left
                wasVisited[index] = true #should be unnecessary now
                index = left
                skip = keys[left].skip
                left = keys[left].left              
                
            elseif (left == 0) && (skip != 0) #is a leaf node, proceed through skip
                #do something
                wasVisited[index] = true
                index = skip
                left = keys[skip].left
                skip = keys[skip].skip
                
                
            elseif (skip != 0) && (left != 0) # is non-right branch, new traversal on skip proceed through left
                wasVisited[index] = true
                index = left
                recursive_traversal(skip, keys, wasVisited, spec)
                skip = keys[left].skip
                left = keys[left].left
                
                
            else
                error("Unexpected traversal exit, datadump:",
                      "index $index, left $left, skip $skip") 
            end
        end
    end
end
    
function is_traversable(keys::Vector{GridKey{T, K}}, spec::SpheresBVHSpecs{T, K}; ShowLonelyKeys=false) where {T, K}
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

function batch_build_traverse(runs::Int, position::Vec3D{T}, spec::SpheresBVHSpecs{T, K}, clct::GenericRandomCollector{T};) where {T, K}
    goodTrees = 0
    selfishTrees = 0
    badTrees = 0

    for each in 1:runs
        treeData = create_mortoncodes(position, spec, clct)
        keys = treeData[1][] 
        I = [GridKey{T, K}(0, 0, MVector{3, K}(0.0, 0.0, 0.0), MVector{3, K}(0.0, 0.0, 0.0), 0, 0) for i in 1:spec.branches_count]
        append!(keys, I)
        update_stackless_bvh!(keys, spec)
        result = is_traversable(keys, spec)
        if result[1]
            goodTrees += 1
        elseif result[2] == "Self-referential"
            selfishTrees += 1
        elseif result[2] == "Nontraversable"
            badTrees += 1
        end
    end

        return goodTrees, selfishTrees, badTrees
 
end
myCollector8 = GenericRandomCollector(; floattype=Float32,
                                    objectnumber=8,
                                    minDim=tuple(0.0, 0.0, 0.0),
                                    maxDim=tuple(1.0, 1.0, 1.0),
                                    temperature=0.01,
                                    randomvelocity=false,
                                    minmass=1.0,
                                    maxmass=5.0,
                                    minimumdistance=0.001,
                                    mincharge=-1f-9,
                                    maxcharge=1f-9
)
position8 =[MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), 
            MVector{3, Float32}(0.346, 0.98, 0.12), MVector{3, Float32}(0.01, 0.76, 0.99), 
            MVector{3, Float32}(0.1111, 0.4, 0.31), MVector{3, Float32}(0.234, 0.29, 0.2), 
            MVector{3, Float32}(0.11346, 0.918, 0.1276), MVector{3, Float32}(0.061, 0.76, 0.989)
]
bvhspec8 = SpheresBVHSpecs(; floattype=Float32, 
                            critical_distance=1.0, 
                            leaves_count=length(position8) 
)
treeData = build_bvh(position8, bvhspec8, myCollector8)
keys = treeData[1][]


##### functionality tests
@testset "traversability" begin


    # "was every leaf (except the last) given a skip rope value?"
    @test keys[8].skip == 0
    for each in 1:7
        @test keys[each].skip != 0
    end

    # "can we repeatedly build traversable trees?"
    goodRuns = batch_build_traverse(100, position8, bvhspec8, myCollector8)
    @test goodRuns[1] == 100 

    # "can tree traversal find every pairing between every atom?"
    list = neighbor_traverse(keys, position8, bvhspec8)
    
end

@testset "bounding volumes" begin
    # "does the root contain the whole scene?"
    for each in 1:8
        for i in eachindex(position8[1])
            @test keys[9].min[i] < position8[each][i] < keys[9].max[i]
        end
    end
end
# @testset "parent boundaries include all atoms" begin
#     for each in 1:7
#         for i in eachindex(position8[1])
#             @test keys[8].min[i] < position8[each][i] < keys[8].max[i]
#         end
#     end
# end

# @testset "morton sorting" begin

# end