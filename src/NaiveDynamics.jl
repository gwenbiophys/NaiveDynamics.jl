module NaiveDynamics

using BenchmarkTools
using UUIDs
using CSV
using DataFrames
using NamedArrays

include("MDInput.jl")
include("Simulator.jl")
include("Testing.jl")
println(GenericRandomCollector(12, 1.0, 2.0, 2.0, 2.1))


"""
    collect_objects(Collector::GenericRandomCollector)

Return a GenericObjectCollection with positions and speeds randomly seeded, as specified by the Collector object

"""

function collect_objects(Collector::GenericRandomCollector)
    positionRange = Uniform(Collector.minsize, Collector.maxsize)
    velocityRange = Uniform(Collector.minspeed, Collector.maxspeed)
    arrayDimensions = (Collector.objectnumber, Collector.objectnumber, Collector.objectnumber)
    
    myObjectCollection = GenericObjectCollection(
        [],
        rand(positionRange, arrayDimensions),
        rand(velocityRange, arrayDimensions),
        fill(uuid1(), Collector.objectnumber)
    )
    return myObjectCollection
end


end # module NaiveDynamics
