using Revise


using BenchmarkTools

using NaiveDynamics
using JLD2
using StaticArrays

using CellListMap


#### no forces simulation
function simulate_noforces(; position, duration, thresh, usebtime=true )
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

    #this was ridiculous to debug. Even though I was closing the jld2 file BEFORE
    #sending it to this funciton, that IO data was still being overwritten.

    bvhspec = SpheresBVHSpecs(;
                                neighbor_distance=thresh, 
                                atom_count=length(position),
                                floattype=Float32
    )
    simspec = SimSpec(; inttype=Int64,
                        floattype=Float32,
                        duration=duration,
                        stepwidth=1,
                        currentstep=1,
                        logLength=10,
                        vDamp=1,
                        threshold=thresh
    )
    #println(clxn.position)

    atoms = clct.objectnumber
    println(" $atoms atoms, $duration duration, $thresh thresh")
    if usebtime
        println("    naive:")
        a = @btime simulate_naive!($clxn, $simspec, $clct)
        println("    bvh:")
        b = @btime simulate_bvh!($clxn, $simspec, $bvhspec, $clct)
    else
        println("    naive:")
        a = @time simulate_naive!(clxn, simspec, clct)
        println("    bvh:")
        b = @time simulate_bvh!(clxn, simspec, bvhspec, clct)
    end
    return nothing
end

f = jldopen("data/positions/positions.jld2", "r")
myposition = deepcopy(read(f, "pos100"))
close(f)

#simulate_noforces(position=myposition, duration=2, thresh=0.03, usebtime=false)
clct = GenericRandomCollector(; floattype=Float32,
                            objectnumber=length(myposition),
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
clxn = collect_objects(clct; position=myposition)

#this was ridiculous to debug. Even though I was closing the jld2 file BEFORE
#sending it to this funciton, that IO data was still being overwritten.
#pos100 breaks at nD 0.15 and atoms perleaf 4
bvhspec = SpheresBVHSpecs(; neighbor_distance=0.1,
                            atom_count=length(myposition),
                            floattype=Float32, 
                            atomsperleaf = 5 
)
simspec = SimSpec(; inttype=Int64,
                        floattype=Float32,
                        duration=10,
                        stepwidth=1,
                        currentstep=1,
                        logLength=10,
                        vDamp=1,
                        threshold=0.03
)




#simulate_noforces(atoms=1024, duration=10, thresh=0.03, usebtime=true)
#Note: These tests are not validated to be consistent, but will hopefulyl show a promising progression through the work
#perm
# 1024 atoms, 10 duration, 0.03 thresh
# 108.138 ms (17683 allocations: 49.85 MiB)
# done wwith naive
# 131.970 ms (21268 allocations: 1.39 MiB)
#direc
# 1024 atoms, 10 duration, 0.03 thresh
# 108.929 ms (17683 allocations: 49.85 MiB)
# done wwith naive
# 148.514 ms (21157 allocations: 1.36 MiB)

#     Jan 8th:
#   naive:
# 80.898 ms (17683 allocations: 49.85 MiB)
#   bvh:
# 87.138 ms (21157 allocations: 1.36 MiB)
#   parallel bvh:
# 32.902 ms (21705 allocations: 1.41 MiB)

#    Jan9th
# 1024 atoms, 10 duration, 0.03 thresh
#   naive:
# 54.374 ms (18145 allocations: 49.89 MiB)
#   bvh:
# 69.047 ms (21223 allocations: 1.70 MiB)
#   parallel bvh:
# 20.867 ms (21761 allocations: 1.75 MiB)




function profile_simulate_noforces(; position, duration, thresh, allocs=false,  allocrate=0.0001)
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
    clxn = collect_objects(clct; position)

    #this was ridiculous to debug. Even though I was closing the jld2 file BEFORE
    #sending it to this funciton, that IO data was still being overwritten. Oh my fault, I didnt copy over, I only read directly

    bvhspec = SpheresBVHSpecs(; neighbor_distance=thresh, 
                                atom_count=length(position),
                                floattype=Float32 
    )
    simspec = SimSpec(; inttype=Int64,
                            floattype=Float32,
                            duration=duration,
                            stepwidth=1,
                            currentstep=1,
                            logLength=10,
                            vDamp=1,
                            threshold=thresh
    )
    #println(clxn.position)

    atoms = clct.objectnumber
    println(" $atoms atoms, $duration duration, $thresh thresh")
    println("    naive:")
    if allocs
        @profview_allocs simulate_naive!(clxn, simspec, clct) sample_rate=allocrate
        println("    bvh:")
        @profview_allocs simulate_bvh!(clxn, simspec, bvhspec, clct) sample_rate=allocrate
    else
        @profview simulate_naive!(clxn, simspec, clct)
        println("    bvh:")
        @profview simulate_bvh!(clxn, simspec, bvhspec, clct)
    end


    return nothing
end

#profile_simulate_noforces(position=myposition, duration=200, thresh=0.03, allocs=false, allocrate=0.001)



function bvh_naive(position, spec, clct; usebtime=true)

    MVecposition = [MVector{3, Float32}(position[i]) for i in eachindex(position)]
    if usebtime
        if length(position) < 10001
            println("  naive:")
            @btime threshold_pairs(unique_pairs($MVecposition), $spec.neighbor_distance)
        end
        a = threshold_pairs(unique_pairs(MVecposition), spec.neighbor_distance)

        println("    bvh:")
        b = @btime build_traverse_bvh($position, $spec)
        
        println("    expt_bvh:")
        c = @btime exptbuild_traverse_bvh($position, $spec)

        println(" CLM.jl:")
        d = @btime neighborlist($position, $spec.neighbor_distance)

        sort!(a, by = x -> x[1])
        sort!(b, by = x -> x[1])
        sort!(c, by = x -> x[1])
        sort!(d, by = x -> x[1])
        println(length(a), " ", length(b), " ", length(c), " ", length(d), " ")

        if length(a) < 55
            println(a)
            println(b)
            println(c)
            println(d)
        end

    else
        # println("    naive:")
        # #a = @time threshold_pairs(unique_pairs(position), spec.neighbor_distance)
        # println("    bvh:")
        # b = build_traverse_bvh(position, spec)
        if length(position) < 10001
            println("  naive:")
            @time threshold_pairs(unique_pairs(MVecposition), spec.neighbor_distance)
        end
        a = threshold_pairs(unique_pairs(MVecposition), spec.neighbor_distance)

        println("    bvh:")
        b = @time build_traverse_bvh(position, spec)
        
        #println("    expt_bvh:")
        #c = @time exptbuild_traverse_bvh(position, spec)

        println(" CLM.jl:")
        d = @time neighborlist(position, spec.neighbor_distance)

        sort!(a, by = x -> x[1])
        sort!(b, by = x -> x[2])
        sort!(b, by = x -> x[1])
        #sort!(c, by = x -> x[1])
        sort!(d, by = x -> x[1])
        #println(length(a), " ", length(b), " ", length(c), " ", length(d), " ")
        println(length(a), " ", length(b), " ", " ", length(d), " ")
        if a != b
            println("BVH and AllToAll are unaligned. Try Again.")
        end
        #println("The ultimate sucess, did I win?: ", a == b)
        if length(a) < 55
            println(a)
            println()
            println(b)

            println()
            #println(c)
            println(d)
        end
    end
end
bvh_naive(myposition, bvhspec, clct; usebtime=false)
# naive:
# 449.664 ms (100 allocations: 897.85 MiB)
#   bvh:
# 139.689 ms (15245 allocations: 10.50 MiB)

# 67.926 μs (0 allocations: 0 bytes)
# 67.436 μs (1 allocation: 16 bytes)
# 29.407 μs (0 allocations: 0 bytes)
# 29.448 μs (0 allocations: 0 bytes)


function bbuild_traverse(position, spec, clct; usebtime=true)
    treeData = TreeData(position, spec)
    #neighbor_traverse(keys, position, spec)
    if usebtime
        println("    base trav:")
        a = @btime neighbor_traverse($treeData.tree, $position, $spec)
        println("    experimental trav:")
        b = @btime expt_neighbor_traverse($treeData.tree, $position, $spec)
    else
        println("    base trav:")
        a = @time neighbor_traverse(keys, position, spec)
        println("    experimental trav:")
        b = @time expt_neighbor_traverse(keys, position, spec)
    end
    return nothing
end

function profile_build_traverse( runs, position, bvhspec, clct; allocs=false, allocrate=0.0001)

    if allocs
        @profview_allocs for i in 1:runs
            build_traverse_bvh(position, bvhspec)
        end
    else
        @profview for i in 1:runs
            build_traverse_bvh(position, bvhspec)
        end
    end
end

#profile_build_traverse(20, myposition, bvhspec, clct)
#From the perspective of the prime thread/the thread which does all of the work
# we spend about 67% performing parallel traversal, 
#about 15% managing the multithreading (most of this time is spend waiting), 
# about 2% sewing the thread work together,
# 3% preparing to build and building the tree,
# about 10% growing the ends of arrays, but i have no idea where that is distributed
# and 3% managing tasks in the Julia runtime.

#Traversal is further split into about 42% of the time performing overlap testing, 12% proximity testing, some memory access outside of those, and
# about 3% that is unaccounted for (points to Line zero of the file)


#atoms=10000, duration=100, thresh=0.03 will fill buffer before completion

#5000 atoms, 100 duration, 0.03 thresh
# run 1:
    #naive: 54% update_pairslist!, 43% threshold_pairs, 2% in unique_pairs, rest in other computations
    #bvh:  53% overlap_test, with an apparent 4% in the function call itself, 8% proximity_test, the  rest in data access from neighbor_traverse
#run 2:
    #naive: 54% update_pairslist!, 44% threshold_pairs, 2% in unique_pairs and the rest
    #bvh: 60% overlap_test, 8% neighbor traverse, the rest in data access and uncertain

#println(myposition[1])
#bbuild_traverse(myposition, bvhspec, clct; usebtime=true)


function exptbuildtraverse(myposition, bvhspec)
    println("   base build")
    #tree = @btime TreeData($myposition, $bvhspec)
    tree =  @time TreeData(myposition, bvhspec)
    println("   expt build")
    #expttree = @btime exptTreeData($myposition, $bvhspec)
    expttree = @time exptTreeData(myposition, bvhspec)

    println("   base trav")
    a = @btime neighbor_traverse($tree.tree, $myposition, $bvhspec)
    println("   expt trav")
    b = @btime expt_neighbor_traverse($expttree.tree, $myposition, $bvhspec)
    freetree = deepcopy(expttree.tree)
    #b = @code_warntype expt_neighbor_traverse(expttree.tree, myposition, bvhspec)
    # println("   expt trav free")
    # b = @btime expt_neighbor_traverse($freetree, $myposition, $bvhspec)
    println(a==b)

    badpos = [MVector{3, Float32}(myposition[i]) for i in eachindex(myposition)]
    c = @btime threshold_pairs(unique_pairs($badpos), $bvhspec.neighbor_distance)

    # d = @btime neighborlist($myposition, $bvhspec.neighbor_distance)
    println(a==c)
    println(b==c)
    println(" base neighbor traverse ", a)
    println(" expt neighbor traverse ", b)
    println(" classic pairlist ", c)
    # println(" CellListMap.jl ", d)



    return nothing
end
#exptbuildtraverse(myposition, bvhspec)
#expttree = @time exptTreeData(myposition, bvhspec)
#b = @code_lowered expt_neighbor_traverse(expttree.tree, myposition, bvhspec)

#treeData = TreeData(myposition, bvhspec)
# keys = treeData.tree
#@benchmark neighbor_traverse($keys, $myposition, $bvhspec)
#@benchmark expt_neighbor_traverse($keys, $myposition, $bvhspec)

#@benchmark build_traverse_bvh($myposition, $bvhspec)


#code native and code llvm look like garbage in the overlap test
#build_traverse_bvh(myposition, bvhspec)






















# println("    naive:")
# a = @btime simulate_naive!($clxn, $simspec, $clct)
# println("    bvh:")
# b = @btime simulate_bvh!($clxn, $simspec, $bvhspec, $clct)
# println("    parallel bvh:")
#b = @benchmark simulate_bvh!($clxn, $simspec, $bvhspec, $clct)