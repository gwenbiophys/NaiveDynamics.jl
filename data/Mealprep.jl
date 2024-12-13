using Revise


using BenchmarkTools
#using CSV
#using StaticArrays
#using NaiveDynamics
using StaticArrays

###### Video recording
using GLMakie

    myCollector = GenericRandomCollector(; floattype=Float32,
                                        objectnumber=30,
                                        minDim=tuple(-1.0, -1.0, -1.0),
                                        maxDim=tuple(1.0, 1.0, 1.0),
                                        temperature=0.01,
                                        randomvelocity=false,
                                        minmass=1.0,
                                        maxmass=5.0,
                                        minimumdistance=0.001,
                                        mincharge=-1f-9,
                                        maxcharge=1f-9
                                        )

    myCollection = collect_objects(myCollector)
    #mySpec = GenericSpec{Int64, Float32}(50, 1, 1, 10, 1)
    mySpec = GenericSpec(; inttype=Int64,
                        floattype=Float32,
                        duration=1000,
                        stepwidth=1,
                        currentstep=1,
                        logLength=10,
                        vDamp=1)
    logpos = simulate!(myCollection, mySpec, myCollector)

    #@profview simulate!(myCollection, mySpec, myCollector)
    #@btime logpos2 = simulate!($myCollection, $mySpec, $myCollector)
    direc = "/home/gwenk/Coding/Julia/NaiveDynamics.jl/data/newhope.mp4"
    record_video(direc, logpos, myCollector; frameinterval = 1)



#using ProfileView #doesnt work in VSCode




#myCollector = GenericRandomCollector{Float32}(40, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -0.000002, 0.000002, 0.001)
#myCollection = collect_objects(myCollector)
#mySpec = GenericSpec{Int64}(800000, 1, 1, 10)
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
#                                     maxcharge=1f-9
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

# mySpec = GenericSpec(; inttype=Int64,
#                     floattype=Float32,
#                     duration=3,
#                     stepwidth=1,
#                     currentstep=1,
#                     logLength=10,
#                     vDamp=1
# )
# @btime simulate!($myCollection, $mySpec, $myCollector1)
# ##### end force testing


###### For BVH
# myCollector2 = GenericRandomCollector(; floattype=Float32,
#                                     objectnumber=3,
#                                     minDim=tuple(0.0, 0.0, 0.0),
#                                     maxDim=tuple(1.0, 1.0, 1.0),
#                                     temperature=0.01,
#                                     randomvelocity=false,
#                                     minmass=1.0,
#                                     maxmass=5.0,
#                                     minimumdistance=0.0001,
#                                     mincharge=-1f-9,
#                                     maxcharge=1f-9
# )
# myCollection1 = collect_objects(myCollector2)

# position8 = [MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), MVector{3, Float32}(0.346, 0.98, 0.12), MVector{3, Float32}(0.01, 0.76, 0.99), MVector{3, Float32}(0.1111, 0.4, 0.31), MVector{3, Float32}(0.234, 0.29, 0.2), MVector{3, Float32}(0.11346, 0.918, 0.1276), MVector{3, Float32}(0.061, 0.76, 0.989) ]
# position7 = [MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), MVector{3, Float32}(0.346, 0.98, 0.12), MVector{3, Float32}(0.01, 0.76, 0.99), MVector{3, Float32}(0.1111, 0.4, 0.31), MVector{3, Float32}(0.234, 0.29, 0.2), MVector{3, Float32}(0.11346, 0.918, 0.1276)]
# position6 = [MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), MVector{3, Float32}(0.346, 0.98, 0.12), MVector{3, Float32}(0.01, 0.76, 0.99), MVector{3, Float32}(0.1111, 0.4, 0.31), MVector{3, Float32}(0.234, 0.29, 0.2)]
# position5 = [MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), MVector{3, Float32}(0.346, 0.98, 0.12), MVector{3, Float32}(0.01, 0.76, 0.99), MVector{3, Float32}(0.1111, 0.4, 0.31)]
# position4 = [MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), MVector{3, Float32}(0.346, 0.98, 0.12), MVector{3, Float32}(0.01, 0.76, 0.99)]
# position3 = [MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), MVector{3, Float32}(0.346, 0.98, 0.12)]

# bvhspec = SpheresBVHSpecs(; floattype=Float32, 
#                             critical_distance=0.3, 
#                             leaves_count=length(position3) 
# )
# # for each in eachindex(position7)
# #     println(minimum(position7[each]))
# # end

# #batch_build_traverse(100, position8, bvhspec, myCollector2, printARun=true)
# #batched_batch_build(10, 100, myCollection1.position, bvhspec, myCollector2)

# build_bvh(position3, bvhspec, myCollector2 )
# bvh_list = build_traverse_bvh(position3, bvhspec, myCollector2)
# naive_list = threshold_pairs(unique_pairs(position3), bvhspec.critical_distance)


##### end bvh


#myCollector = GenericRandomCollector(Float32, 40, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -0.02, 0.02, 0.001)
#myCollection = collect_objects(myCollector)
#mySpec = GenericSpec{Int64}(4000, 1, 1, 10)
#logpos2 = simulate_dumloop!(myCollection, mySpec, myCollector)
#@btime logpos2 = simulate!($myCollection, $mySpec, $myCollector)
