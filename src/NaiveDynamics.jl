module NaiveDynamics

using BenchmarkTools
using Revise

include("MDInput.jl")
include("Simulator.jl")

# These should all be made into tests
"""
testCollector = GenericCollector(3, 0.0, 10.0, -2.0, 2.0)
println("Success for testCollector! ", testCollector)
myTestCollection = collect_objects(testCollector)
println("Here are your positions: ", myTestCollection.position)
println("Position dimensions number: ", ndims(myTestCollection.position))
#println(length(myTestCollection.uniqueID))

mySystem = GenericSystem(10, 1, 0, myTestCollection)
mySimulation = GenericSimulation(mySystem, true)
simulate!(mySimulation)
println("Here are your positions after the simulation: ", myTestCollection.position)
println("Position dimensions after number: ", ndims(myTestCollection.position))
"""
end # module NaiveDynamics
