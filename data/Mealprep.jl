using Revise


#using BenchmarkTools
#using CSV
#using StaticArrays
using NaiveDynamics
using StaticArrays
import GLMakie

    # myCollector = GenericRandomCollector(; floattype=Float32,
    #                                     objectnumber=50,
    #                                     minDim=tuple(-1.0, -1.0, -1.0),
    #                                     maxDim=tuple(1.0, 1.0, 1.0),
    #                                     temperature=0.01,
    #                                     randomvelocity=false,
    #                                     minmass=1.0,
    #                                     maxmass=5.0,
    #                                     minimumdistance=0.001,
    #                                     mincharge=-1f-9,
    #                                     maxcharge=1f-9
    #                                     )

    # myCollection = collect_objects(myCollector)
    # #mySpec = GenericSpec{Int64, Float32}(50, 1, 1, 10, 1)
    # mySpec = GenericSpec(; inttype=Int64,
    #                     floattype=Float32,
    #                     duration=100,
    #                     stepwidth=1,
    #                     currentstep=1,
    #                     logLength=10,
    #                     vDamp=1)
    # logpos = simulate!(myCollection, mySpec, myCollector)
    # #@profview simulate!(myCollection, mySpec, myCollector)
    # #@btime logpos2 = simulate!($myCollection, $mySpec, $myCollector)
    # direc = "/home/gwenk/Coding/Julia/NaiveDynamics.jl/data/iWant.mp4"
    # record_video(direc, logpos, myCollector; frameinterval = 1)



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



# myCollector1 = GenericRandomCollector(; floattype=Float32,
#                                     objectnumber=8,
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
myCollector2 = GenericRandomCollector(; floattype=Float32,
                                    objectnumber=4,
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
myCollection1 = collect_objects(myCollector2)

position8 = [MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), MVector{3, Float32}(0.346, 0.98, 0.12), MVector{3, Float32}(0.01, 0.76, 0.99), MVector{3, Float32}(0.1111, 0.4, 0.31), MVector{3, Float32}(0.234, 0.29, 0.2), MVector{3, Float32}(0.11346, 0.918, 0.1276), MVector{3, Float32}(0.061, 0.76, 0.989) ]
position7 = [MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), MVector{3, Float32}(0.346, 0.98, 0.12), MVector{3, Float32}(0.01, 0.76, 0.99), MVector{3, Float32}(0.1111, 0.4, 0.31), MVector{3, Float32}(0.234, 0.29, 0.2), MVector{3, Float32}(0.11346, 0.918, 0.1276)]
position6 = [MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), MVector{3, Float32}(0.346, 0.98, 0.12), MVector{3, Float32}(0.01, 0.76, 0.99), MVector{3, Float32}(0.1111, 0.4, 0.31), MVector{3, Float32}(0.234, 0.29, 0.2)]
position5 = [MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), MVector{3, Float32}(0.346, 0.98, 0.12), MVector{3, Float32}(0.01, 0.76, 0.99), MVector{3, Float32}(0.1111, 0.4, 0.31)]
position4 = [MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), MVector{3, Float32}(0.346, 0.98, 0.12), MVector{3, Float32}(0.01, 0.76, 0.99)]
position3 = [MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), MVector{3, Float32}(0.346, 0.98, 0.12)]

bvhspec = SpheresBVHSpecs(; floattype=Float32, 
                            interaction_distance=0.1, 
                            leaves_count=length(position4) 
)

function wellwell(a, b)
    countsa = 0
    countsb = 0
    for i in 1:100000
        if a < b
            countsa += 1
        else
            countsb += 1
        end
    end
    println(countsa, "<-a ", countsb)
end
#@code_native wellwell(5, 5)
#batch_build_traverse(100, position4, bvhspec, myCollector2, printARun=true)

build_bvh(position4, bvhspec, myCollector2 )
# function dist(pos1, pos2)
#     d = 0.0
#     for e in eachindex(pos2)
#         d1 = (pos1[e] - pos2[e])^2
#         d += d1
#     end
#     println(sqrt(d))
# end
# dist(position[1], position[2])
# dist(position[2], position[3])

#myCollection1.

#myCollector = GenericRandomCollector(Float32, 40, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -0.02, 0.02, 0.001)
#myCollection = collect_objects(myCollector)
#mySpec = GenericSpec{Int64}(4000, 1, 1, 10)
#logpos2 = simulate_dumloop!(myCollection, mySpec, myCollector)
#@btime logpos2 = simulate!($myCollection, $mySpec, $myCollector)
