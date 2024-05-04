#module MDInput
# File for generating inputs to the "simulation!" functions
# and eventually input-file reading from custom and/or community-based tooling


# 1. Generic input-generation, last updated 4/08/2024




export 
# would also have to include definitions of types and functions
    ObjectCollection,
    GenericObjectCollection,
    Collector,
    GenericRandomCollector,
    GenericZeroCollector,
    GenericUserValueCollector,
    collect_objects,
    testCollector,
    myTestCollection



abstract type Collector end
"""
    GenericRandomCollector(objectnumber, minsize, maxsize, minspeed, maxspeed)
    
An Collector-subtype meant to acquire additional information from the user about how to make their system.
In the collection function, positions and velocities will be randomly seeded from this Collector's boundary values.
    
"""



abstract type ObjectCollection end
mutable struct GenericObjectCollection <: ObjectCollection 
    currentstep::AbstractArray{Integer,1}
    name::AbstractArray{String, 1}
    #mass::AbstractArray{Number, 1}
    #radius::AbstractArray{AbstractFloat, 1}
    index::AbstractArray{Integer, 1}
    position::AbstractArray{SizedVector{3, AbstractFloat}, 1}
    velocity::AbstractArray{SizedVector{3, AbstractFloat}, 1}
    force::AbstractArray{SizedVector{3, AbstractFloat}, 1}

    #uniqueID::AbstractArray{UUID,1}

end




"""
    Collector

Collector super-type for simulation initialization.

"""
mutable struct GenericRandomCollector <: Collector
    objectnumber::Integer
    minsize::AbstractFloat
    maxsize::AbstractFloat
    minspeed::AbstractFloat
    maxspeed::AbstractFloat
end
mutable struct GenericZeroCollector <: Collector
    objectnumber::Integer
end
mutable struct GenericUserValueCollector <: Collector
    objectnumber::Integer
    userposition::AbstractFloat
    uservelocity::AbstractFloat
end

"""
    collect_objects(Collector::GenericRandomCollector)

Return a GenericObjectCollection with positions and speeds randomly seeded, as specified by the Collector object

"""

function collect_objects(Collector::GenericRandomCollector)
    positionRange = Uniform(Collector.minsize, Collector.maxsize)
    velocityRange = Uniform(Collector.minspeed, Collector.maxspeed)

    objectcount = Collector.objectnumber
    step_n=1

    simCollection = GenericObjectCollection(
        fill(step_n, objectcount),
        fill("duck", objectcount),
        [1:objectcount;],
        [SizedVector{3, Float64}(rand(positionRange, 3)) for each in 1:objectcount],
        [SizedVector{3, Float64}(rand(velocityRange, 3)) for each in 1:objectcount],
        [SizedVector{3, Float64}(zeros(Float64, 3)) for each in 1:objectcount],
        )
    return simCollection
end

function collect_objects(Collector::GenericZeroCollector)
    arrayDimensions = (Collector.objectnumber, Collector.objectnumber, Collector.objectnumber)

    objectnumber = Collector.objectnumber
    step_n=1
    myObjectCollection = DataFrame(
        current_step=step_n,
        object_index= collect(Int, 1:objectnumber),
        position_x=zeros(Float64, length(1:objectnumber)),
        position_y=zeros(Float64, length(1:objectnumber)),
        position_z=zeros(Float64, length(1:objectnumber)),
        velocity_x=zeros(Float64, length(1:objectnumber)),
        velocity_y=zeros(Float64, length(1:objectnumber)),
        velocity_z=zeros(Float64, length(1:objectnumber))
        )
    return myObjectCollection
end

function collect_objects(Collector::GenericUserValueCollector )
    #arrayDimensions = (Collector.objectnumber, Collector.objectnumber, Collector.objectnumber)
    position = Collector.userposition
    velocity = Collector.uservelocity
    objectnumber = Collector.objectnumber
    step_n=1
    myObjectCollection = DataFrame(
        current_step=step_n,
        object_index= collect(Int, 1:objectnumber),
        position_x=fill(position, length(1:objectnumber)),
        position_y=fill(position, length(1:objectnumber)),
        position_z=fill(position, length(1:objectnumber)),
        velocity_x=fill(velocity, length(1:objectnumber)),
        velocity_y=fill(velocity, length(1:objectnumber)),
        velocity_z=fill(velocity, length(1:objectnumber))
        )

    return myObjectCollection
end

#end #module
