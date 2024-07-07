using Revise
using BenchmarkTools
using CSV
using StaticArrays
using GLMakie

#using ProfileView #doesnt work in VSCode




myCollector = GenericRandomCollector(Float32, 40, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -0.00002, 0.00002, 0.001)
myCollection = collect_objects(myCollector)
mySpec = GenericSpec(40000, 1, 1, 10)
#@btime simulate!($largeSimulation, $largeCollector)
log = simulate_bravado!(myCollection, mySpec, myCollector)

#println("log length ", length(log))
#println(typeof(log))

record_video("./NaiveDynamics.jl/data/newTest.mp4", log, myCollector )


function hope()
    a = GenericRandomCollector(Float32, 100, -0.5, -0.6, -0.7, 0.5, 0.6, 0.7, -0.002, 0.002, 0.001)
    b::GenericStaticRandomCollector = GenericStaticRandomCollector(Float32, 100, -0.5, -0.6, -0.7, 0.5, 0.6, 0.7, -0.002, 0.002, 0.001)
    CollectorVec = [a]

    println("Mutable Vecs!")
    largeCollector = a
    largeCollection = collect_objects(largeCollector)
    largeSystem::GenericSystem = GenericSystem(2000, 1, 1, largeCollection)
    largeSimulation::GenericSimulation = GenericSimulation(largeSystem, 5)

    println("   Generator expressions")
    @btime simulate_MVec!($largeSimulation, $largeCollector)
    println("   Loop exclusion")
    @btime simulate!($largeSimulation, $largeCollector)
    #@code_warntype simulate!(largeSimulation, largeCollector)

    println("Static Vecs!")
    sCollector = b
    sCollection = collect_objects(sCollector)
    sSystem::GenericSystem = GenericSystem(2000, 1, 1, sCollection)
    sSimulation::GenericSimulation = GenericSimulation(sSystem, 5)

    println("   Generator expressions")
    @btime simulate_SVec!($sSimulation, $sCollector)
    println("   Loop exclusion")
    #@code_warntype simulate!($largeSimulation, $largeCollector)
    @btime simulate!($sSimulation, $sCollector)
    println()

end
#hope()





#largeCollector = GenericRandomCollector(Float32, 1000, -0.5, -0.6, -0.7, 0.5, 0.6, 0.7, -0.002, 0.002, 0.001)
#largeCollection = collect_objects(largeCollector)
#largeSystem::GenericSystem = GenericSystem(2000000, 1, 1, largeCollection)
#largeSimulation::GenericSimulation = GenericSimulation(largeSystem, 5)
#@profview simulate_MVec!(largeSimulation, largeCollector)

#println(largeHistory.Simulation.collection.position)
# 28.200 μs (3057 allocations: 122.11 KiB) with mutable Vector
# 13.600 μs (57 allocations: 78.73 KiB) static vector

#mutable genericObjColl: 4.3 ms, 706.78 KiB, 19631 allocs
#immutable:                 4.23 ms, 707.41 KiB, 19639 allocs
# really interesting!

#function batchrun(a)
#    result = nothing
#    for i in  1:1
#        Collector = GenericRandomCollector(Float32, 10, -0.05, -0.06, -0.07, 0.05, 0.06, 0.07, -0.2, 0.2)
#        Collection = collect_objects(Collector)
#        System = GenericSystem(3, 1, 1, Collection)
#        simulation = GenericSimulation(System, true)
#        
#        result = simulate!(simulation, Collector)
#  
#
#        
#    end
#    return result
#end
#batchrun(c)

#println(simlog)
#CSV.write("large_rand.csv", largeCollection)
# with dataFrames and 100 000 particles, simulation fails due to memory error
