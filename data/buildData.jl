using Revise
using BenchmarkTools
using CSV
using StaticArrays
using GLMakie

#using ProfileView #doesnt work in VSCode




myCollector = GenericRandomCollector(Float32, 40, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -0.002, 0.002, 0.001)
myCollection = collect_objects(myCollector)
mySpec = GenericSpec{Int64}(4000, 1, 1, 10)
#@btime simulate!($largeSimulation, $largeCollector)
#@btime logpos = simulate_bravado!($myCollection, $mySpec, $myCollector)
#@code_warntype simulate_bravado!(myCollection, mySpec, myCollector)
logpos = simulate!(myCollection, mySpec, myCollector)
#@btime logpos2 = simulate!($myCollection, $mySpec, $myCollector)

record_video("./NaiveDynamics.jl/data/newhope.mp4", logpos, myCollector; frameinterval = 10 )


#myCollector = GenericRandomCollector(Float32, 40, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -0.02, 0.02, 0.001)
#myCollection = collect_objects(myCollector)
#mySpec = GenericSpec{Int64}(4000, 1, 1, 10)
#logpos2 = simulate_dumloop!(myCollection, mySpec, myCollector)
#@btime logpos2 = simulate!($myCollection, $mySpec, $myCollector)
