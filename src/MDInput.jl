#module MDInput
# File for generating inputs to the "simulation!" functions
# and eventually input-file reading from custom and/or community-based tooling


# 1. Generic input-generation, last updated 4/08/2024




export 
# would also have to include definitions of types and functions
    Vec3D,
    StatVec3D,
    ObjectCollection,
    GenericObjectCollection,
    GenericStaticCollection,
    Collector,
    GenericRandomCollector,
    GenericStaticRandomCollector,
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

const Vec3D{T} = Vector{MVector{3, T}}
const StatVec3D{T} = Vector{SVector{3, T}}

abstract type ObjectCollection end
mutable struct GenericObjectCollection{K, L, T<:AbstractFloat} <: ObjectCollection 
    currentstep::Vector{K}
    name::Vector{String}
    mass::Vector{L}
    radius::Vector{T}
    index::Vector{K}
    position::Vector{MVector{3, T}}
    velocity::Vector{MVector{3, T}}
    force::Vector{MVector{3, T}}

    #uniqueID::Vector{UUID,1}

end

struct GenericStaticCollection{K, L, T<:AbstractFloat} <: ObjectCollection 
    currentstep::Vector{K}
    name::Vector{String}
    mass::Vector{L}
    radius::Vector{T}
    index::Vector{K}
    position::Vector{SVector{3, T}}
    velocity::Vector{SVector{3, T}}
    force::Vector{SVector{3, T}}

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
    
    
    # i am uncertain if the 3D vector or the 1D vector will be easier.
    #minDimension::MVector{1, AbstractFloat}
    #maxDimension::MVector{1, AbstractFloat}
    #minDimension::MVector{3, AbstractFloat}
    #maxDimension::MVector{3, AbstractFloat}

    minspeed::T
    maxspeed::T
    minimumdistance::T
end

struct GenericStaticRandomCollector{T<:AbstractFloat} <: Collector
    T
    objectnumber::Integer
    min_xDim::T
    min_yDim::T
    min_zDim::T
    max_xDim::T
    max_yDim::T
    max_zDim::T
    
    
    # i am uncertain if the 3D vector or the 1D vector will be easier. by i will go 1D
    #minDimension::MVector{1, AbstractFloat}
    #maxDimension::MVector{1, AbstractFloat}
    #minDimension::MVector{3, AbstractFloat}
    #maxDimension::MVector{3, AbstractFloat}

    minspeed::T
    maxspeed::T
    minimumdistance::T
end
#function generate_randomcollector(objectcount, min_xDim, min_yDim, min_zDim, max_xDim, max_yDim, max_zDim, minspeed, maxspeed)
    #return GenericRandomCollector(objectcount, SA[min_xDim, min_yDim, min_zDim], SA[max_xDim, max_yDim, max_zDim,], minspeed, maxspeed)
#end

struct GenericZeroCollector <: Collector
    objectnumber::Integer
end

struct GenericUserValueCollector{T<:AbstractFloat} <: Collector

    userposition::AbstractArray{MVector{3, T}, 1}
    uservelocity::AbstractArray{MVector{3, T}, 1}

    T
    objectnumber::Integer
    min_xDim::T
    min_yDim::T
    min_zDim::T
    max_xDim::T
    max_yDim::T
    max_zDim::T
    
    
    # i am uncertain if the 3D vector or the 1D vector will be easier. by i will go 1D
    #minDimension::MVector{1, AbstractFloat}
    #maxDimension::MVector{1, AbstractFloat}
    #minDimension::MVector{3, AbstractFloat}
    #maxDimension::MVector{3, AbstractFloat}

    minspeed::T
    maxspeed::T
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

    xyz = [MVector{3, Float32}(x[i], y[i], z[i]) for i in 1:objectcount]

    return xyz
end
function generate_positions(Collector::GenericStaticRandomCollector)
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

    xyz = [SVector{3, Float32}(x[i], y[i], z[i]) for i in 1:objectcount]

    return xyz
end
function generate_onePosition(Collector::GenericStaticRandomCollector)
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
    return SVector{3, Float32}(x, y, z)
end

"""
    function generate_pruned_positions(Collector::GenericRandomCollector, Collection)

After collect_objects() has generated a vector of vectors of positions, 
this function will naively prune the positions and replace them with other positions. 
At present, pruning proceeds until try-catch control flow until neighborlist()
fails due to having no neighbors within the minimiumdistance cutoff.
"""

function unique_pairs_prune(a::AbstractArray, threshold::AbstractFloat)
    # TODO only push unique pairs to the list for eachindex(a), instead of for each pair
    tooClose = 0
    counter = 0
    j_cutoff = length(a)-1

    
    dx = 1.0
    dy = 1.0
    dz = 1.0
    d2 = 1.0


    #does not even theoretically work until this line is deleted. just needs 1 more gloss over with the brain
    for i in 1:length(a)-1
            for j in i+1:length(a)-1 
                dx = a[i][1] - a[j][1]
                dy = a[i][2] - a[j][2] 
                dz = a[i][3] - a[j][3]
                d2 = sqrt(dx^2 + dy^2 + dz^2)  
                if d2 < threshold
                    fill!(a[i], Inf64)
                    tooClose += 1

                end   
            end

    end
    return tooClose
end

function generate_pruned_positions(Collector::GenericRandomCollector, Collection)
    #1. generate a set of radialPositions, in each each ooooh no wait i'd have to make radii spheres
    #2. ask if the distance between any of the spheres is negative with a neighborlist cutoff =o (if that works)
    minDist = Collector.minimumdistance
    
    #try
        #presumably this implementation will proceed indefinitely until it crashes into an error.
        #Hopefully the computer doesn't spend too long on whether this is an error
        
        tooClose = unique_pairs_prune(Collection.position, minDist;)
        while tooClose < 0
            # to be deleted once i see it's working
            #for each in eachindex(tooClose)
                #we want to create a function that will zero-out/Inf-out any position set 
                #which has an index sitting in the neighborlist
                #we want to replace the position values labeled by index
                #for and in eachindex(Collection.position)
                    #if tooClose[and][1] == Collection.position[and]
                        #for every in eachindex(Collection.position[and])
                            #every = Inf64 # TODO change to match the precision selected by user
                        #end
                    #end
                #end
            #end
            for object in eachindex(Collection.position)
                if Collection.position[object] == Inf64 # TODO change to NaN ?
                    
                    fill!(Collection.position[object], generate_onePosition(Collector))
                end
            end
            tooClose = unique_pairs_prune(Collection.position, minDist;)
        end
        return # i dont think this is necessary, as this is a mutating function with no output, but just in case
            
    #catch
        #return
    #end

    
end

struct SimDomainLacksDensity <: Exception
    message::String
end
function density_check(collection::GenericObjectCollection, collector::GenericRandomCollector)
    # TODO turn this into a special error and/or get rid of this function once neighborlist() does not fail
    LJCutoff_neighborList = neighborlist(collection.position, 0.5)
    if length(LJCutoff_neighborList) == 0
       println(throw("The simulation domain is not dense enough for CellListMap to work, please add more particles and/or reduce box size. 
       Otherwise, the simulation may only proceed through each step until neighborlist() throws an OutOfMemory error."))
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
    simCollection = GenericObjectCollection{Int64, Float64, Collector.T}(
        fill(step_n, objectcount),
        fill("duck", objectcount),
        rand(1:5, objectcount),
        fill(0.01, objectcount),
        [1:objectcount;],
        generate_positions(Collector),
        [MVector{3, Collector.T}(rand(velocityRange, 3)) for each in 1:objectcount],
        [MVector{3, Collector.T}(zeros(Float64, 3)) for each in 1:objectcount],
        )
    #density_check(simCollection, Collector)

    generate_pruned_positions(Collector, simCollection)
    return simCollection
end

function collect_objects(Collector::GenericStaticRandomCollector)

    velocityRange = Uniform(Collector.minspeed, Collector.maxspeed)

    objectcount = Collector.objectnumber
    step_n=1

    # TODO update the Collector Series so that the user can input names and masses and radii.
    simCollection = GenericStaticCollection(
        fill(step_n, objectcount),
        fill("duck", objectcount),
        rand(1:5, objectcount),
        fill(0.01, objectcount),
        [1:objectcount;],
        generate_positions(Collector),
        [SVector{3, Collector.T}(rand(velocityRange, 3)) for each in 1:objectcount],
        [SVector{3, Collector.T}(zeros(Float32, 3)) for each in 1:objectcount],
        )
    #these methods are nonfunctional with SVectors
    #density_check(simCollection, Collector)
    #generate_pruned_positions(Collector, simCollection)
    return simCollection
end

function collect_objects(Collector::GenericZeroCollector)
    arrayDimensions = (Collector.objectnumber, Collector.objectnumber, Collector.objectnumber)

    objectnumber = Collector.objectnumber
    step_n=1
    myObjectCollection = DataFrame(
        current_step=step_n,
        object_index= collect(Int, 1:objectnumber),
        position_x=zeros(Float32, length(1:objectnumber)),
        position_y=zeros(Float32, length(1:objectnumber)),
        position_z=zeros(Float32, length(1:objectnumber)),
        velocity_x=zeros(Float32, length(1:objectnumber)),
        velocity_y=zeros(Float32, length(1:objectnumber)),
        velocity_z=zeros(Float32, length(1:objectnumber))
        )
    
    return myObjectCollection
end

function collect_objects(Collector::GenericUserValueCollector )
    #arrayDimensions = (Collector.objectnumber, Collector.objectnumber, Collector.objectnumber)
    position = Collector.userposition
    velocity = Collector.uservelocity
    objectcount = Collector.objectnumber
    step_n=1
    myObjectCollection = GenericObjectCollection(
        fill(step_n, objectcount),
        fill("duck", objectcount),
        rand(1:5, objectcount),
        fill(0.01, objectcount),
        [1:objectcount;],
        position,
        velocity,
        [MVector{3, Collector.T}(zeros(Float32, 3)) for each in 1:objectcount],
        )

    return myObjectCollection
end

#end #module
