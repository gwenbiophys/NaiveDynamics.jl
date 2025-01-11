using Revise


using BenchmarkTools

using NaiveDynamics
using JLD2


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
                                        maxcharge=1f-9
    )
    clxn = collect_objects(clct)

    #this was ridiculous to debug. Even though I was closing the jld2 file BEFORE
    #sending it to this funciton, that IO data was still being overwritten.
    copyto!.(clxn.position, position)
    bvhspec = SpheresBVHSpecs(; floattype=Float32, 
                                critical_distance=thresh, 
                                leaves_count=length(position) 
    )
    simspec = GenericSpec(; inttype=Int64,
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
myposition = deepcopy(read(f, "pos5000"))
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
                            maxcharge=1f-9
)
clxn = collect_objects(clct)

#this was ridiculous to debug. Even though I was closing the jld2 file BEFORE
#sending it to this funciton, that IO data was still being overwritten.
copyto!.(clxn.position, myposition)
bvhspec = SpheresBVHSpecs(; floattype=Float32, 
                            critical_distance=0.2, 
                            leaves_count=length(myposition) 
)
simspec = GenericSpec(; inttype=Int64,
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
    maxcharge=1f-9
)
    clxn = collect_objects(clct)

    #this was ridiculous to debug. Even though I was closing the jld2 file BEFORE
    #sending it to this funciton, that IO data was still being overwritten. Oh my fault, I didnt copy over, I only read directly
    copyto!.(clxn.position, position)
    bvhspec = SpheresBVHSpecs(; floattype=Float32, 
                                critical_distance=thresh, 
                                leaves_count=length(position) 
    )
    simspec = GenericSpec(; inttype=Int64,
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





function bcreate_mortoncodes(; atoms, thresh=0.01)
    clct = GenericRandomCollector(; floattype=Float32,
                                        objectnumber=atoms,
                                        minDim=tuple(0.0, 0.0, 0.0),
                                        maxDim=tuple(1.0, 1.0, 1.0),
                                        temperature=0.01,
                                        randomvelocity=false,
                                        minmass=1.0,
                                        maxmass=5.0,
                                        minimumdistance=0.0001,
                                        mincharge=-1f-9,
                                        maxcharge=1f-9
    )
    clxn = collect_objects(clct)
    bvhspec = SpheresBVHSpecs(; floattype=Float32, 
                                critical_distance=thresh, 
                                leaves_count=length(clxn.position) 
    )
    println(atoms, " atoms")
    println("permsort")
    @time create_mortoncodes(clxn.position, bvhspec, clct)
    println("direct sort")
    @time create_mortoncodes_direct(clxn.position, bvhspec, clct)
    println("now btime")
    println("permsort")
    @btime create_mortoncodes($clxn.position, $bvhspec, $clct)
    println("direct sort")
    @btime create_mortoncodes_direct($clxn.position, $bvhspec, $clct)
    return nothing
end

function bbuild_traverse(position, spec, clct; usebtime=true)
    treeData = build_bvh(position, spec, clct)
    keys = treeData[1][]
    #neighbor_traverse(keys, position, spec)
    if usebtime
        println("    base:")
        a = @btime neighbor_traverse($keys, $position, $spec)
        println("    experimental:")
        b = @btime expt_neighbor_traverse($keys, $position, $spec)
    else
        println("    base:")
        a = @time neighbor_traverse(keys, position, spec)
        println("    experimental:")
        b = @time expt_neighbor_traverse(keys, position, spec)
    end
    return nothing
end




#atoms=10000, duration=100, thresh=0.03 will fill buffer before completion

#5000 atoms, 100 duration, 0.03 thresh
# run 1:
    #naive: 54% update_pairslist!, 43% threshold_pairs, 2% in unique_pairs, rest in other computations
    #bvh:  53% overlap_test, with an apparent 4% in the function call itself, 8% proximity_test, the  rest in data access from neighbor_traverse
#run 2:
    #naive: 54% update_pairslist!, 44% threshold_pairs, 2% in unique_pairs and the rest
    #bvh: 60% overlap_test, 8% neighbor traverse, the rest in data access and uncertain
#bcreate_mortoncodes(atoms=1024)

#println(myposition[1])
bbuild_traverse(myposition, bvhspec, clct; usebtime=true)

# treeData = build_bvh(myposition, bvhspec, clct)
# keys = treeData[1][]
#@benchmark neighbor_traverse($keys, $myposition, $bvhspec)
#@benchmark expt_neighbor_traverse($keys, $myposition, $bvhspec)

































# println("    naive:")
# a = @btime simulate_naive!($clxn, $simspec, $clct)
# println("    bvh:")
# b = @btime simulate_bvh!($clxn, $simspec, $bvhspec, $clct)
# println("    parallel bvh:")
#b = @benchmark simulate_bvh!($clxn, $simspec, $bvhspec, $clct)