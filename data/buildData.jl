module buildData
using BenchmarkTools
using CSV
include("C:/Users/Trist/.julia/dev/NaiveMD/NaiveDynamics.jl/src/NaiveDynamics.jl")
include("C:/Users/Trist/.julia/dev/NaiveMD/NaiveDynamics.jl/src/MDInput.jl")
include("C:/Users/Trist/.julia/dev/NaiveMD/NaiveDynamics.jl/src/Simulator.jl")

# 1. Build a  simulation with user-locked values.
valuesCollector = GenericUserValueCollector(10, 0.2, 12.9)
valuesCollection = collect_objects(valuesCollector)
valuesSystem = GenericSystem(5, 1, 1, valuesCollection)
valuesSimulation = GenericSimulation(valuesSystem, true)
theResult = simulate!(valuesSimulation)
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
largeCollector = GenericRandomCollector(100000, -5.0, 5.0, -0.2, 0.2)
largeCollection = collect_objects(largeCollector)
largeSystem = GenericSystem(1, 1, 1, largeCollection)
largeSimulation = GenericSimulation(largeSystem, true)
@btime simulate!(largeSimulation)
#println(largeCollection)
#CSV.write("large_rand.csv", largeCollection)
# with dataFrames and 100 000 particles, simulation fails due to memory error
end