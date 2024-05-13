using Revise
using BenchmarkTools
using CSV



#include("C:/Users/Trist/.julia/dev/NaiveMD/NaiveDynamics.jl/src/NaiveDynamics.jl")
#include("C:/Users/Trist/.julia/dev/NaiveMD/NaiveDynamics.jl/src/MDInput.jl")
#include("C:/Users/Trist/.julia/dev/NaiveMD/NaiveDynamics.jl/src/Simulator.jl")

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

#largeCollector = GenericRandomCollector(Float64, 102, -5.0, -6.0, -7.0, 5.0, 6.0, 7.0, -0.2, 0.2)
#largeCollection = collect_objects(largeCollector)
#largeSystem = GenericSystem(10, 1, 1, largeCollection)
#largeSimulation = GenericSimulation(largeSystem, true)
#largeHistory = simulate!(largeSimulation, largeCollector)
#println(largeHistory.position)


function batchrun(a)
    result = nothing
    for i in  1:200
        Collector = GenericRandomCollector(Float64, i, -5.0, -6.0, -7.0, 5.0, 6.0, 7.0, -0.2, 0.2)
        Collection = collect_objects(Collector)
        System = GenericSystem(10, 1, 1, Collection)
        simulation = GenericSimulation(System, true)
        
        result = simulate!(simulation, Collector)
  

        
    end
    return result
end
batchrun(c)

#println(simlog)
#CSV.write("large_rand.csv", largeCollection)
# with dataFrames and 100 000 particles, simulation fails due to memory error
