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
    #generate_randomcollector,
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
mutable struct GenericObjectCollection{K, L, T<:AbstractFloat} <: ObjectCollection 
    currentstep::AbstractArray{K,1}
    name::AbstractArray{String, 1}
    mass::AbstractArray{L, 1}
    radius::AbstractArray{T, 1}
    index::AbstractArray{K, 1}
    position::AbstractArray{MVector{3, T}, 1}
    velocity::AbstractArray{MVector{3, T}, 1}
    force::AbstractArray{MVector{3, T}, 1}

    #uniqueID::AbstractArray{UUID,1}

end




"""
    Collector

Collector super-type for simulation initialization.

"""
struct GenericRandomCollector{T<:AbstractFloat} <: Collector
    T
    objectnumber::Integer
    min_xDim::T
    min_yDim::T
    min_zDim::T
    max_xDim::T
    max_yDim::T
    max_zDim::T
    
    
    # i am uncertain if the 3D vector or the 1D vector will be easier. by i will go 1D
    #minDimension::SVector{1, AbstractFloat}
    #maxDimension::SVector{1, AbstractFloat}
    #minDimension::SVector{3, AbstractFloat}
    #maxDimension::SVector{3, AbstractFloat}

    minspeed::T
    maxspeed::T
end
#function generate_randomcollector(objectcount, min_xDim, min_yDim, min_zDim, max_xDim, max_yDim, max_zDim, minspeed, maxspeed)
    #return GenericRandomCollector(objectcount, SA[min_xDim, min_yDim, min_zDim], SA[max_xDim, max_yDim, max_zDim,], minspeed, maxspeed)
#end

struct GenericZeroCollector <: Collector
    objectnumber::Integer
end

struct GenericUserValueCollector <: Collector
    objectnumber::Integer
    userposition::AbstractFloat
    uservelocity::AbstractFloat
end

"""
    generate_positions(objectcount, min_xDim, min_yDim, min_zDim, max_xDim, max_yDim, max_zDim)

Return a position vector of 3-dimensional mutable vectors when given box dimensions from  the GenericRandomCollector object.
"""

function generate_positions(Collector::GenericRandomCollector)
    # TODO maybe jsut close these in. for development i built a function outside of this file and then
    # fitted the function here, calling on the stack of equal signs. but a quick copy paste could make it briefer
    min_xDim = Collector.min_xDim
    min_yDim = Collector.min_yDim
    min_zDim = Collector.min_zDim
    max_xDim = Collector.max_xDim
    max_yDim = Collector.max_yDim
    max_zDim = Collector.max_zDim
    objectcount = Collector.objectnumber

    xDimRange = Uniform(min_xDim, max_xDim)
    yDimRange = Uniform(min_yDim, max_yDim)
    zDimRange = Uniform(min_zDim, max_zDim)

    x = rand(xDimRange, objectcount)
    y = rand(yDimRange, objectcount)
    z = rand(zDimRange, objectcount)

    xyz = [MVector{3, Float64}(x[i], y[i], z[i]) for i in 1:objectcount]

    return xyz
end
function generate_prune_positions(Collector::GenericRandomCollector, position, radius)
    #1. generate a set of radialPositions, in each each ooooh no wait i'd have to make radii spheres
    #2. ask if the distance between any of the spheres is negative with a neighborlist cutoff =o (if that works)
    
    radialPosition = copy(position)
    for i in eachindex(position)
        radialPosition[i] .+= radius[i]
        for i in eachindex(position)
            if radialPosition[i]
            
            end

        end
    end
    
end


"""
    collect_objects(Collector::GenericRandomCollector)

Return a GenericObjectCollection with positions and speeds randomly seeded, as specified by the Collector object

"""

function collect_objects(Collector::GenericRandomCollector)

    velocityRange = Uniform(Collector.minspeed, Collector.maxspeed)

    objectcount = Collector.objectnumber
    step_n=1

    # TODO update the Collector Series so that the user can input names and masses and radii.
    simCollection = GenericObjectCollection(
        fill(step_n, objectcount),
        fill("duck", objectcount),
        rand(1:5, objectcount),
        fill(0.01, objectcount),
        [1:objectcount;],
        generate_positions(Collector),
        [MVector{3, Collector.T}(rand(velocityRange, 3)) for each in 1:objectcount],
        [MVector{3, Collector.T}(zeros(Float64, 3)) for each in 1:objectcount],
        )
    #generate_prune_positions(Collector, simCollection.position, simCollection.radius)
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
