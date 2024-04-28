#module MDInput
# File for generating inputs to the "simulation!" functions
# and eventually input-file reading from custom and/or community-based tooling


# 1. Generic input-generation, last updated 4/08/2024

using UUIDs
using Distributions
using DataFrames


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

abstract type ObjectCollection end
mutable struct GenericObjectCollection <: ObjectCollection 
    name::AbstractArray{String, 1}
    #mass::AbstractArray{Number, 1}
    #radius::AbstractArray{AbstractFloat, 1}
    position::AbstractArray{AbstractFloat, 3}
    velocity::AbstractArray{AbstractFloat, 3}
    #totalobjects
    uniqueID::AbstractArray{UUID,1}
    #forcevectors
end

# this in an incorrect method i am certain, as there is no AbstractNamedArray object, 
# but I want a struct to look like this rather than just instantiating a global variable to hold the specific type
# though maybe i could just instantiate this NamedArray at a function call, like how other objects are currently
# instantiated.
#mutable struct GenericObjectCollection2 <: NamedArray
 #   name::AbstractArray{String, 1}
#    mass::AbstractArray{Number, 1}
 #   radius::AbstractArray{AbstractFloat, 1}
 #   position::AbstractArray{AbstractFloat, 3}
 #   velocity::AbstractArray{AbstractFloat, 3}
  #  uniqueID::AbstractArray{UUID,1}
#end

#GenericObjectCollection3 = NamedArray(
   # AbstractArray{String, 1},
  #  AbstractArray{Number, 1},
  #  AbstractArray{AbstractFloat, 1},
   # AbstractArray{AbstractFloat, 3},
   # AbstractArray{AbstractFloat, 3},
   # AbstractArray{AbstractFloat, 3},
   # AbstractArray{UUID,1};
 #   "name",
 #   "mass",
  #  "radius",
  #  "position",
   # "velocity",
   # "force",
   # "uniqueID"
 #   )



"""
    Collector

Collector super-type for simulation initialization.

"""
abstract type Collector end
"""
    GenericRandomCollector(objectnumber, minsize, maxsize, minspeed, maxspeed)

An Collector-subtype meant to acquire additional information from the user about how to make their system.
In the collection function, positions and velocities will be randomly seeded from this Collector's boundary values.

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
    arrayDimensions = (Collector.objectnumber, Collector.objectnumber, Collector.objectnumber)
    
    myObjectCollection = GenericObjectCollection(
        [],
        rand(positionRange, arrayDimensions),
        rand(velocityRange, arrayDimensions),
        fill(uuid1(), Collector.objectnumber)
    )
    return myObjectCollection
end

# is there a more julian way of doing this?
function collect_objects(Collector::GenericZeroCollector)
    arrayDimensions = (Collector.objectnumber, Collector.objectnumber, Collector.objectnumber)
    
    myObjectCollection = GenericObjectCollection(
        [],
        zeros(Float64, arrayDimensions),
        zeros(Float64, arrayDimensions),
        fill(uuid1(), Collector.objectnumber)
    )
    return myObjectCollection
end

function collect_objects(Collector::GenericUserValueCollector, uservalue::Float64)
    arrayDimensions = (Collector.objectnumber, Collector.objectnumber, Collector.objectnumber)
    
    myObjectCollection = GenericObjectCollection(
        [],
        fill(uservalue, arrayDimensions),
        fill(uservalue, arrayDimensions),
        fill(uuid1(), Collector.objectnumber)
    )
    return myObjectCollection
end

#end #module
