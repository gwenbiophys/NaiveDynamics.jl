using Revise
using BenchmarkTools
using CSV
using StaticArrays
using Cthulhu #doesnt work in VSCode

#### Testing for the 4 bugs when running neighborlist()
# 1. Test for directly overlapping particles

#positions = [MVector{3, Float64}(fill(1.0, 3)) for each in 1:5]
#velocities = [MVector{3, Float64}(zeros(3)) for each in 1:5]
#Collector = GenericUserValueCollector(positions, velocities, Float64, 5, -5.0, -6.0, -7.0, 5.0, 6.0, 7.0, -0.2, 0.2)
#Collection = collect_objects(Collector)
#System = GenericSystem(100, 1, 1, Collection)
#simulation = GenericSimulation(System, true)
#result = simulate!(simulation, Collector)
## truncation error, forces go to infinty, positions go to infinity, nobody's happy


#### Benchtesting and eventually generating starting data
# 1. Build a  simulation with user-locked values.
#valuesCollector = GenericUserValuesCollector(10, 0.2, 12.9)
#valuesCollection = collect_objects(valuesCollector)
#valuesSystem = GenericSystem(5, 1, 1, valuesCollection)
#valuesSimulation = GenericSimulation(valuesSystem, true)
#theResult = simulate!(valuesSimulation)
#println("And here is the result!: ", theResult)
#resultantTable = DataFrame(step=theResult.steparray, positions=theResult.positionrecord, velocities=theResult.velocityrecord)
#println(theResult)

# 2. Small, random GenericDataFrame
"""
smallCollector = GenericRandomCollector(10, -5.0, 5.0, -0.2, 0.2)
smallCollection = collect_objects(smallCollector)
CSV.write("small_rand.csv", smallCollection)
"""

# 3. ENHANCE
a = [25, 100, 1000, 10000]
b = [1, 2, 3, 4, 5, 10, 15, 20, 25, 100]
c = [500]
function hope()
    a = GenericRandomCollector(Float64, 100, -0.5, -0.6, -0.7, 0.5, 0.6, 0.7, -0.002, 0.002, 0.001)
    b = GenericStaticRandomCollector(Float64, 100, -0.5, -0.6, -0.7, 0.5, 0.6, 0.7, -0.002, 0.002, 0.001)
    CollectorVec = [a]
    for i in eachindex(CollectorVec)
        largeCollector = CollectorVec[i]
        largeCollection = collect_objects(largeCollector)
        largeSystem = GenericSystem(10, 1, 1, largeCollection)
        largeSimulation = GenericSimulation(largeSystem, 5)
        @btime simulate!($largeSimulation, $largeCollector)
        #if i == 1
            #@btime simulate_unified!($largeSimulation, $largeCollector)
       # else
            
        #end
        #@btime simulate_unified!($largeSimulation, $largeCollector)
        #@btime simulate_oneloop!($largeSimulation, $largeCollector)
    end
end
hope()
#println(largeHistory.Simulation.collection.position)
# 28.200 μs (3057 allocations: 122.11 KiB) with mutable Vector
# 13.600 μs (57 allocations: 78.73 KiB) static vector

#mutable genericObjColl: 4.3 ms, 706.78 KiB, 19631 allocs
#immutable:                 4.23 ms, 707.41 KiB, 19639 allocs
# really interesting!

#function batchrun(a)
#    result = nothing
#    for i in  1:1
#        Collector = GenericRandomCollector(Float64, 10, -0.05, -0.06, -0.07, 0.05, 0.06, 0.07, -0.2, 0.2)
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
