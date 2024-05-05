using Revise
using BenchmarkTools
using CSV
using NaiveDynamics


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
largeCollector = GenericRandomCollector(10000, -5.0, 5.0, -0.2, 0.2)
largeCollection = collect_objects(largeCollector)
largeSystem = GenericSystem(1, 1, 1, largeCollection)
largeSimulation = GenericSimulation(largeSystem, true)

@btime simulate!(largeSimulation, largeCollector)
# 1 step at 100 000  particles with push!= 109.952ms, 2596967 and 56.41 MiB
    # with no call to record_simulation_bench = 50.535 ms (800003 allocations: 15.26 MiB)
# 1 step at 100k with append! = 
#println(simlog)
#CSV.write("large_rand.csv", largeCollection)
# with dataFrames and 100 000 particles, simulation fails due to memory error
