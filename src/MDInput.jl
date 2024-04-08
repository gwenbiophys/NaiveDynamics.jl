#module MDInput
# File for generating inputs to the "simulation!" functions
# and eventually input-file reading from custom and/or community-based tooling


# 1. Generic input-generation, last updated 3/22/2024
using Revise
using UUIDs
using Distributions


export 
# would also have to include definitions of types and functions
    ObjectCollection,
    GenericObjectCollection,
    Collector,
    GenericCollector,
    collect_objects,
    testCollector,
    myTestCollection


abstract type ObjectCollection end
mutable struct GenericObjectCollection <: ObjectCollection 
    name::AbstractArray{String, 1}
    #mass::AbstractArray{Number, 1}
    #radius::AbstractArray{AbstractFloat, 1}
    position::AbstractArray{AbstractFloat, 3}
    velocity::AbstractArray{AbstractFloat, 3}
    #totalobjects
    uniqueID::AbstractArray{UUID,1}
end

abstract type Collector end
mutable struct GenericCollector <: Collector
    objectnumber::Integer
    minsize::AbstractFloat
    maxsize::AbstractFloat
    minspeed::AbstractFloat
    maxspeed::AbstractFloat
end

function collect_objects(Collector::GenericCollector)
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






function genericFunction()
    return println("GoodJob!")
end



function genericVectorFill!(vectorName,data)
    fill!(vectorName, data)

end

#end #module
