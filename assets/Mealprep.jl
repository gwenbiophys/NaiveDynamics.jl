using Revise


using BenchmarkTools
#using CSV
#using StaticArrays
using NaiveDynamics
using StaticArrays
using JET

###### Video recording
using GLMakie

    myCollector = GenericRandomCollector(; floattype=Float32,
                                        objectnumber=6,
                                        minDim=tuple(-1.0, -1.0, -1.0),
                                        maxDim=tuple(1.0, 1.0, 1.0),
                                        temperature=0.01,
                                        randomvelocity=false,
                                        minmass=1.0,
                                        maxmass=5.0,
                                        minimumdistance=0.001,
                                        mincharge=-1f-9,
                                        maxcharge=1f-9,
                                        pregeneratedposition=false
                                        )

    myCollection = collect_objects(myCollector)
    #mySpec = SimSpec{Int64, Float32}(50, 1, 1, 10, 1)
    mySpec = SimSpec(; inttype=Int64,
                        floattype=Float32,
                        duration=500,
                        stepwidth=1,
                        currentstep=1,
                        logLength=10,
                        vDamp=1,
                        threshold = 0.1
    )
    logpos = simulate!(myCollection, mySpec, myCollector)

    #@profview simulate!(myCollection, mySpec, myCollector)
    #@btime logpos2 = simulate!($myCollection, $mySpec, $myCollector)
    direc = "data/newhope.mp4"#"/home/gwenk/Coding/Julia/NaiveDynamics.jl/data/newhope.mp4"
    record_video(direc, logpos, myCollector; frameinterval = 1)



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

bvhspec = SpheresBVHSpecs(; bounding_distance=0.003, 
                            neighbor_distance=0.003, 
                            leaves_count=length(myCollection1.position),
                            floattype=Float32 
)
simspec = SimSpec(; inttype=Int64,
                    floattype=Float32,
                    duration=60,
                    stepwidth=1,
                    currentstep=1,
                    logLength=10,
                    vDamp=1,
                    threshold=0.03
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
    treeData = TreeData(position, bvhspec)
    #keys = treeData[1][]
    for i in 1:runs

        neighbor_traverse(treeData[1][], position, bvhspec)
        TreeData!(treeData, position, bvhspec)

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
#safe_position = [IndexSafePosition{Float32, Int32}(i, SVector{3, Float32}(position[i])) for i in eachindex(position)]
mutsafe_position = [MutableIndexSafePosition{Float32, Int32}(i, position[i]) for i in eachindex(position)]
indy = [0 for each in eachindex(position)]

# @btime sortperm!($indy, $position, by=x->[1])
# @btime sort!($position, by=x->x[1])
# @btime sort!($sposition, by=x -> x[1])
# @btime sort!($safe_position, by=x -> x.vec[1])
# @btime sort!($mutsafe_position, by=x -> x.vec[1])
#Base.Threads.lock

#TreeData(position, bvhspec)
#build_traverse_bvh(position, bvhspec)
#@benchmark threshold_pairs(unique_pairs($position), $bvhspec.neighbor_distance)
#list = @btime build_traverse_bvh($position, $bvhspec)
#naivelist = @btime threshold_pairs(unique_pairs($position), $bvhspec.neighbor_distance)

#g = run_bvh(1, position, bvhspec)
# at 10 000 leaves, only run <10 times
# at ~<1000 leaves, run 1000 times
# h = @btime   run_bvh(10, $position, $bvhspec) # at 100 objects: 271.214 ms (239976 allocations: 23.92 MiB)
# d = @btime run_naive(10, $position, $bvhspec.neighbor__distance)

# i = @profview   run_bvh(40000, position, bvhspec, myCollector2)
# j = @profview_allocs   run_bvh(40000, position, bvhspec, myCollector2) sample_rate = 0.001 #default 0.0001

#  g = @btime   run_bvh(1000, $position, $bvhspec, $myCollector2) # at 100 objects: 271.214 ms (239976 allocations: 23.92 MiB)
#  d = @btime run_naive(1000, $position, $bvhspec.neighbor_distance)
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
# bvhspec8 = SpheresBVHSpecs(; bounding_distance=10.0, 
#                             neighbor_distance=10.0, 
#                             leaves_count=length(position8),
#                             floattype=Float32,
# )
#TreeData(position8, bvhspec8)
#@btime with vectorized bounding volume update = 2.002 μs (123 allocations: 5.31 KiB)
#Btime with forloop bounding volume updating = 1.954 μs (116 allocations: 5.09 KiB)

#bvh_list = @btime build_traverse_bvh($position, $bvhspec)
#naive_list = @btime threshold_pairs(unique_pairs($position), $bvhspec.neighbor_distance)


#@time create_mortoncodes(position, bvhspec, myCollector2)
#@time create_mortoncodes_perm(position, bvhspec, myCollector2)
#@profview TreeData(position, bvhspec)
#@profview_allocs TreeData(position, bvhspec)
#@btime TreeData($position, $bvhspec)
#bvh_list = build_traverse_bvh(position, bvhspec) 

#run_naive(10, position, bvhspec.neighbor_distance)
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





######## SIMD overlap_test experiments
# Base.@propagate_inbounds function newoverlap_test(each, k, min, pos, max)

#     counter = 0
    
#     for i in eachindex(pos[k])
#         yes=0
#         low=false
#         high=false

#         if min[each][i] < pos[k][i]
#             yes+=1
#             #low=true
#         end

#         if pos[k][i] < max[each][i]
#             yes+=1
#             #high=true
#         end

#         if yes==2
#             counter +=1
#         end
#         # if low & high
#         #     counter +=1
#         # end

#     end
#     return counter
# end
# Base.@propagate_inbounds function overlap_test(i, k, min, pos, max)
#     #@assert length(min[i]) == length(max[i]) == length(pos[k])
#     return all(min[i] .< pos[k] .< max[i])

# end

# # min = [SVector{3, Float32}(rand(3)) for i in 1:10000]
# # pos = [SVector{3, Float32}(rand(3)) for i in 1:10000]
# # maxi = [SVector{3, Float32}(rand(3)) for i in 1:10000]
# min = [NTuple{3, Float32}(rand(3)) for i in 1:10000]
# pos = [NTuple{3, Float32}(rand(3)) for i in 1:10000]
# maxi = [NTuple{3, Float32}(rand(3)) for i in 1:10000]
# # maxisimd = [Vec{3, Float32}((maxi[i])) for i in 1:10000]
# # minisimd = [Vec{3, Float32}((min[i])) for i in 1:10000]
# # posisimd = [Vec{3, Float32}((pos[i])) for i in 1:10000]
# maxisimd = [Vec{4, Float32}((maxi[i][1], maxi[i][2], maxi[i][3], maxi[i][3])) for i in 1:10000]
# minisimd = [Vec{4, Float32}((min[i][1], min[i][2], min[i][3], min[i][3])) for i in 1:10000]
# posisimd = [Vec{4, Float32}((pos[i][1], pos[i][2], pos[i][3], pos[i][3])) for i in 1:10000]
# i = 3
# k = 5
# subject = Vec{8, Float32}((min[i][1], min[i][2], min[i][3], pos[k][1], pos[k][2], pos[k][3],  Float32(0.0), Float32(0.0)))
# query = Vec{8, Float32}((pos[k][1], pos[k][2], pos[k][3], maxi[i][1], maxi[i][2], maxi[i][3],  Float32(1.0), Float32(1.0)))
# bigsubject = [Vec{8, Float32}((min[i][1], min[i][2], min[i][3], pos[i][1], pos[i][2], pos[i][3],  Float32(0.0), Float32(0.0))) for i in eachindex(pos)]
# bigquery = [Vec{8, Float32}((pos[i][1], pos[i][2], pos[i][3], maxi[i][1], maxi[i][2], maxi[i][3],  Float32(1.0), Float32(1.0))) for i in eachindex(pos)]

# Base.@propagate_inbounds function simdoverlap_interior(subject::Vec{N, T}, query::Vec{N, T}) where {N, T}
#     return all(subject < query)
# end

# Base.@propagate_inbounds function resimdoverlap_interior(subject::Vector{Vec{N, T}}, query::Vector{Vec{N, T}}) where {N, T}
#     sum = 0
#     for each in eachindex(subject) 
#         sum += all(subject[each] < query[each])
#     end
#     return sum
# end


# Base.@propagate_inbounds function simdoverlap_test(i::Int, k::Int, min::Vector{Vec{N, T}}, pos::Vector{Vec{N, T}}, maxi::Vector{Vec{N, T}}) where {N,T}

#     #subject = Vec{8, Float32}((min[i][1], min[i][2], min[i][3], pos[k][1], pos[k][2], pos[k][3], Float32(1.0), Float32(1.0)))
#     #query = Vec{8, Float32}((pos[k][1], pos[k][2], pos[k][3], maxi[i][1], maxi[i][2], maxi[i][3], Float32(0.0), Float32(0.0)))
#     # query = Vec{3, Float32}((pos[k][1], pos[k][2], pos[k][3]))#, pos[k][3])) #duplicated
#     # low = Vec{3, Float32}((min[i][1], min[i][2], min[i][3]))#, min[i][3]))
#     # high = Vec{3, Float32}((maxi[i][1], maxi[i][2], maxi[i][3]))#, maxi[i][3]))

#     # query = Vec{3, Float32}(pos[k])#, pos[k][3])) #duplicated
#     # low = Vec{3, Float32}(min[i])#, min[i][3]))
#     # high = Vec{3, Float32}(maxi[i])#, maxi[i][3]))
#     # query = pos[k]
#     # low = min[i]
#     # high = maxi[i]
#     # c = low < query

#     # d = query < high

#     c = min[i] < pos[k]
#     d = pos[k] < maxi[i]

#     #e = c + d
#     e = c & d

#     #f = Vec{3, Int64}(e)
#     return all(e)
#     #return all(e)
#     #return reduce(&, e)
#     #return sum(NTuple{3, Bool}(e))
#     #return sum(NTuple{3, Bool}(e))
#     #return sum(f)
#     #sum(low < query < high)

# end
# #a = @code_native @inbounds overlap_test(3, 5, min, pos, maxi)
# # f = @btime @inbounds overlap_test(3, 5, $min, $pos, $maxi)
# # a = @btime @inbounds newoverlap_test(3, 5, $min, $pos, $maxi)
# # x = @btime @inbounds sum($min[3] .< $pos[5] .< $maxi[3])
# # b = @btime @inbounds simdoverlap_test(3, 5, $min, $pos, $maxi)
# # x = @time @inbounds sum(min[3] .< pos[5] .< maxi[3])
# # b =  simdoverlap_test(3, 5, min, pos, maxi) 
# Base.@propagate_inbounds function reoverlap(min, pos, maxi, iters)
#     #println(length(pos))
#     sum = 0
#     #each = 1
#     for i in 1:1:iters
#         for each in eachindex(pos)
#         #k = length(pos) - i + 1
#             sum += @inbounds overlap_test(each, each, min, pos, maxi)
#         end
#     end
#     return sum
# end

# Base.@propagate_inbounds function resimd(min, pos, maxi, iters)
#     sum = 0
#     #each = 1
#     for i in 1:1:iters
#         for each in eachindex(pos)
#         k = length(pos) - i + 1
#             sum += @inbounds simdoverlap_test(each, each, min, pos, maxi)

#         end
#     end
#     return sum
# end

# @btime @inbounds overlap_test(3, 5, $min, $pos, $maxi)
# @btime @inbounds simdoverlap_test(3, 5, $minisimd, $posisimd, $maxisimd)
# @btime @inbounds simdoverlap_interior($subject, $query)
# println()
# b = @btime @inbounds reoverlap($min, $pos, $maxi, 1)
# a = @btime @inbounds resimd($minisimd, $posisimd, $maxisimd, 1)
# yesssss = @btime @inbounds resimdoverlap_interior($bigsubject, $bigquery)
# # println(pos[1])
# # println(min[1])
# # println(maxi[1])
#  println(reoverlap(min, pos, maxi, 1))
#  println(resimd(minisimd, posisimd, maxisimd, 1))
#  println(resimdoverlap_interior(bigsubject, bigquery))



# #a = @code_native @inbounds sum(min[3] .< pos[5] .< maxi[3])
# #a = @code_native @inbounds overlap_test(3, 5, min, pos, maxi)
# #println(f," ", a, " ", b)
# i = 3
# k = 5
# subject = Vec{8, Float32}((min[i][1], min[i][2], min[i][3], pos[k][1], pos[k][2], pos[k][3], Float32(1.0), Float32(1.0)))
# query = Vec{8, Float32}((pos[k][1], pos[k][2], pos[k][3], maxi[i][1], maxi[i][2], maxi[i][3], Float32(0.0), Float32(0.0)))
# #b = @code_native @inbounds sum(subject < query)
# println()
# println()
# println()
# #x = @code_native @inbounds sum(min[3] .< pos[5] .< maxi[3])