using Revise


using BenchmarkTools
#using CSV
#using StaticArrays
using NaiveDynamics
using StaticArrays
using JET

###### Video recording
# using GLMakie

#     myCollector = GenericRandomCollector(; floattype=Float32,
#                                         objectnumber=6,
#                                         minDim=tuple(-1.0, -1.0, -1.0),
#                                         maxDim=tuple(1.0, 1.0, 1.0),
#                                         temperature=0.01,
#                                         randomvelocity=false,
#                                         minmass=1.0,
#                                         maxmass=5.0,
#                                         minimumdistance=0.001,
#                                         mincharge=-1f-9,
#                                         maxcharge=1f-9,
#                                         pregeneratedposition=false
#                                         )

#     myCollection = collect_objects(myCollector)
#     #mySpec = SimSpec{Int64, Float32}(50, 1, 1, 10, 1)
#     mySpec = SimSpec(; inttype=Int64,
#                         floattype=Float32,
#                         duration=1000,
#                         stepwidth=1,
#                         currentstep=1,
#                         logLength=10,
#                         vDamp=1)
#     logpos = simulate!(myCollection, mySpec, myCollector)

#     #@profview simulate!(myCollection, mySpec, myCollector)
#     #@btime logpos2 = simulate!($myCollection, $mySpec, $myCollector)
#     direc = "/home/gwenk/Coding/Julia/NaiveDynamics.jl/data/newhope.mp4"
#     record_video(direc, logpos, myCollector; frameinterval = 1)



#using ProfileView #doesnt work in VSCode




#myCollector = GenericRandomCollector{Float32}(40, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -0.000002, 0.000002, 0.001)
#myCollection = collect_objects(myCollector)
#mySpec = SimSpec{Int64}(800000, 1, 1, 10)
# memory limit starts getting feisty at 4 million steps
#@btime simulate!($largeSimulation, $largeCollector)
#@btime logpos = simulate_bravado!($myCollection, $mySpec, $myCollector)
#@code_warntype simulate_bravado!(myCollection, mySpec, myCollector)
#@profview logpos = simulate!(myCollection, mySpec, myCollector)
#myCollector = GenericRandomCollector{Float32}(50, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 0.0001, false, 1.0, 5.0, 0.001, -5.0, 5.0)


# ##### force testing
# myCollector1 = GenericRandomCollector(; floattype=Float32,
#                                     objectnumber=80,
#                                     minDim=tuple(-1.0, -1.0, -1.0),
#                                     maxDim=tuple(1.0, 1.0, 1.0),
#                                     temperature=0.01,
#                                     randomvelocity=false,
#                                     minmass=1.0,
#                                     maxmass=5.0,
#                                     minimumdistance=0.001,
#                                     mincharge=-1f-9,
#                                     maxcharge=1f-9,
#                                      pregeneratedposition = false
# )
# myCollection = collect_objects(myCollector1)

# function force_testing(runs)
#     myCollector1 = GenericRandomCollector(; floattype=Float32,
#                                         objectnumber=800,
#                                         minDim=tuple(-1.0, -1.0, -1.0),
#                                         maxDim=tuple(1.0, 1.0, 1.0),
#                                         temperature=0.01,
#                                         randomvelocity=false,
#                                         minmass=1.0,
#                                         maxmass=5.0,
#                                         minimumdistance=0.001,
#                                         mincharge=-1f-9,
#                                         maxcharge=1f-9
                                        #   pregeneratedposition=false
#     )
#     myCollection = collect_objects(myCollector1)
#     force_C = deepcopy(myCollection.force)
#     pairslist = unique_pairs(myCollection.position)
#     for i in 1:runs
#         @btime force_coulomb!($force_C, $pairslist, $myCollection.charge)
#         force_C = deepcopy(myCollection.force)
#     end
# end
# #force_testing(3)

# mySpec = SimSpec(; inttype=Int64,
#                     floattype=Float32,
#                     duration=3,
#                     stepwidth=1,
#                     currentstep=1,
#                     logLength=10,
#                     vDamp=1
# )
# @btime simulate!($myCollection, $mySpec, $myCollector1)
                                                                        # ##### end force testing


###### Neighbor Search: random days
myCollector2 = GenericRandomCollector(; floattype=Float32,
                                    objectnumber=5000,
                                    minDim=tuple(0.0, 0.0, 0.0),
                                    maxDim=tuple(1.0, 1.0, 1.0),
                                    temperature=0.01,
                                    randomvelocity=false,
                                    minmass=1.0,
                                    maxmass=5.0,
                                    minimumdistance=0.0001,
                                    mincharge=-1f-9,
                                    maxcharge=1f-9,
                                    pregeneratedposition=false
)
myCollection1 = collect_objects(myCollector2)
position = generate_positions(myCollector2)
#println(sizeof(position))


position8 = [MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), MVector{3, Float32}(0.346, 0.98, 0.12), MVector{3, Float32}(0.01, 0.76, 0.99), MVector{3, Float32}(0.1111, 0.4, 0.31), MVector{3, Float32}(0.234, 0.29, 0.2), MVector{3, Float32}(0.11346, 0.918, 0.1276), MVector{3, Float32}(0.061, 0.76, 0.989) ]
position7 = [MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), MVector{3, Float32}(0.346, 0.98, 0.12), MVector{3, Float32}(0.01, 0.76, 0.99), MVector{3, Float32}(0.1111, 0.4, 0.31), MVector{3, Float32}(0.234, 0.29, 0.2), MVector{3, Float32}(0.11346, 0.918, 0.1276)]
position6 = [MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), MVector{3, Float32}(0.346, 0.98, 0.12), MVector{3, Float32}(0.01, 0.76, 0.99), MVector{3, Float32}(0.1111, 0.4, 0.31), MVector{3, Float32}(0.234, 0.29, 0.2)]
position5 = [MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), MVector{3, Float32}(0.346, 0.98, 0.12), MVector{3, Float32}(0.01, 0.76, 0.99), MVector{3, Float32}(0.1111, 0.4, 0.31)]
position4 = [MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), MVector{3, Float32}(0.346, 0.98, 0.12), MVector{3, Float32}(0.01, 0.76, 0.99)]
position3 = [MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), MVector{3, Float32}(0.346, 0.98, 0.12)]

bvhspec = SpheresBVHSpecs(; floattype=Float32, 
                            critical_distance=0.003, 
                            leaves_count=length(myCollection1.position) 
)
simspec = SimSpec(; inttype=Int64,
                    floattype=Float32,
                    duration=60,
                    stepwidth=1,
                    currentstep=1,
                    logLength=10,
                    vDamp=1
)
##### Compare sim time runs
#a = simulate_naive!(myCollection1, simspec, myCollector2)
# println("done wwith naive")
#b = @time simulate_bvh!(myCollection1, simspec, bvhspec, myCollector2)
# for each in eachindex(position7)
#     println(minimum(position7[each]))
# end

#batch_build_traverse(100, position8, bvhspec, myCollector2, printARun=true)
#batched_batch_build(10, 100, myCollection1.position, bvhspec, myCollector2)
function run_bvh(runs, position, bvhspec, myCollector2)
    treeData = build_bvh(position, bvhspec)
    #keys = treeData[1][]
    for i in 1:runs

        neighbor_traverse(treeData[1][], position, bvhspec)
        rebuild_bvh!(treeData, position, bvhspec)

    end

end
function run_naive(runs, position, thresh)
    list = unique_pairs(position)

    for i in 1:runs
        update_pairslist!(position, list)
        threshold_pairs(list, thresh)
    end
end
sort
sposition = [SVector{3, Float32}(position[i]) for i in eachindex(position)]
safe_position = [IndexSafePosition{Float32, Int32}(i, SVector{3, Float32}(position[i])) for i in eachindex(position)]
mutsafe_position = [MutableIndexSafePosition{Float32, Int32}(i, position[i]) for i in eachindex(position)]
indy = [0 for each in eachindex(position)]

# @btime sortperm!($indy, $position, by=x->[1])
# @btime sort!($position, by=x->x[1])
# @btime sort!($sposition, by=x -> x[1])
# @btime sort!($safe_position, by=x -> x.vec[1])
# @btime sort!($mutsafe_position, by=x -> x.vec[1])
#Base.Threads.lock

#build_bvh(position, bvhspec)
build_traverse_bvh(position, bvhspec)
#@benchmark threshold_pairs(unique_pairs($position), $bvhspec.critical_distance)
#list = @btime build_traverse_bvh($position, $bvhspec)
#naivelist = @btime threshold_pairs(unique_pairs($position), $bvhspec.critical_distance)

#g = run_bvh(1, position, bvhspec)
# at 10 000 leaves, only run <10 times
# at ~<1000 leaves, run 1000 times
# h = @btime   run_bvh(10, $position, $bvhspec) # at 100 objects: 271.214 ms (239976 allocations: 23.92 MiB)
# d = @btime run_naive(10, $position, $bvhspec.critical_distance)

# i = @profview   run_bvh(40000, position, bvhspec, myCollector2)
# j = @profview_allocs   run_bvh(40000, position, bvhspec, myCollector2) sample_rate = 0.001 #default 0.0001

#  g = @btime   run_bvh(1000, $position, $bvhspec, $myCollector2) # at 100 objects: 271.214 ms (239976 allocations: 23.92 MiB)
#  d = @btime run_naive(1000, $position, $bvhspec.critical_distance)
#g = run_bvh(1, position, bvhspec, myCollector2)



###### Neighbor Search: Using fixed position
# myCollector8 = GenericRandomCollector(; floattype=Float32,
#                                     objectnumber=8,
#                                     minDim=tuple(0.0, 0.0, 0.0),
#                                     maxDim=tuple(1.0, 1.0, 1.0),
#                                     temperature=0.01,
#                                     randomvelocity=false,
#                                     minmass=1.0,
#                                     maxmass=5.0,
#                                     minimumdistance=0.001,
#                                     mincharge=-1f-9,
#                                     maxcharge=1f-9,
                                    # pregeneratedposition=false
# )
# position8 =[MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), 
#             MVector{3, Float32}(0.346, 0.98, 0.12), MVector{3, Float32}(0.01, 0.76, 0.99), 
#             MVector{3, Float32}(0.1111, 0.4, 0.31), MVector{3, Float32}(0.234, 0.29, 0.2), 
#             MVector{3, Float32}(0.11346, 0.918, 0.1276), MVector{3, Float32}(0.061, 0.76, 0.989)
# ]
# bvhspec8 = SpheresBVHSpecs(; floattype=Float32, 
#                             critical_distance=10.0, 
#                             leaves_count=length(position8) 
# )
#build_bvh(position8, bvhspec8)
#@btime with vectorized bounding volume update = 2.002 μs (123 allocations: 5.31 KiB)
#Btime with forloop bounding volume updating = 1.954 μs (116 allocations: 5.09 KiB)

#bvh_list = @btime build_traverse_bvh($position, $bvhspec)
#naive_list = @btime threshold_pairs(unique_pairs($position), $bvhspec.critical_distance)


#@time create_mortoncodes(position, bvhspec, myCollector2)
#@time create_mortoncodes_perm(position, bvhspec, myCollector2)
#@profview build_bvh(position, bvhspec)
#@profview_allocs build_bvh(position, bvhspec)
#@btime build_bvh($position, $bvhspec)
#bvh_list = build_traverse_bvh(position, bvhspec) 

#run_naive(10, position, bvhspec.critical_distance)
#run_perms(10, position, bvhspec, myCollector2)
#naive_list = 
#println(bvh_list[1])
#println(naive_list[1])

##### end bvh


#myCollector = GenericRandomCollector(Float32, 40, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -0.02, 0.02, 0.001)
#myCollection = collect_objects(myCollector)
#mySpec = SimSpec{Int64}(4000, 1, 1, 10)
#logpos2 = simulate_dumloop!(myCollection, mySpec, myCollector)
#@btime logpos2 = simulate!($myCollection, $mySpec, $myCollector)
