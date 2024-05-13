#module Simulator

#println("Starting up!")
export 
    System,
    Simulation,
    TimeMethod,
    GenericIntegrator,
    Logger,
    GenericLogger,
    GenericSystem,
    GenericSimulation,
    simulate!,
    record_simulation



abstract type System end
abstract type Simulation end
abstract type TimeMethod end
struct GenericIntegrator <: TimeMethod
    stepwidth::Integer
end


abstract type Logger end
mutable struct GenericLogge <: Logger 
    steparray::AbstractArray{Integer, 1}
    positionrecord::AbstractArray{AbstractArray{AbstractFloat,3}, 1}
    velocityrecord::AbstractArray{AbstractArray{AbstractFloat,3}, 1}
end
mutable struct GenericLogger <: Logger 
    currentstep::AbstractArray{AbstractArray{Integer, 1},1}
    name::AbstractArray{AbstractArray{String, 1}, 1}
    #mass::AbstractArray{AbstractArray{Number, 1}, 1}
    #radius::AbstractArray{AbstractFloat, 1}
    index::AbstractArray{AbstractArray{Integer, 1}, 1}
    position::AbstractArray{AbstractArray{MVector{3, AbstractFloat}, 1}, 1}
    velocity::AbstractArray{AbstractArray{MVector{3, AbstractFloat}, 1}, 1}
    force::AbstractArray{AbstractArray{MVector{3, AbstractFloat}, 1}, 1}

    #uniqueID::AbstractArray{UUID,1}

end

mutable struct GenericLoggerAppend <: ObjectCollection 
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

mutable struct PositionLogger <: ObjectCollection 
    currentstep::AbstractArray{Integer,1}
    name::AbstractArray{String, 1}
    index::AbstractArray{Integer, 1}
    position::AbstractArray{SizedVector{3, AbstractFloat}, 1}
end

struct GenericSystem <: System
    duration::Integer
    stepwidth::Integer
    currentstep::Integer
    objectcollection
end
struct GenericSimulation <: Simulation
    system::GenericSystem
    do_logging::Bool
end

#function record_simulation(step_n, position, velocity, myLog)
    
 #   push!(myLog.steparray, step_n)
 #   push!(myLog.positionrecord, position)
   # push!(myLog.velocityrecord, velocity)
  #  return #println("You are logging!")
#end

function record_simulation(simCollection, simlog)
    
    #push!(simCollection, simlog)
    #this is ultimately how the function should work, just a clean shove
    #but GenericObjectCollection may have to become a single array and not a mutable datatype
    #and then we can just push the entire array to the log
    push!(simlog.currentstep, simCollection.currentstep)
    push!(simlog.name, simCollection.name)
    push!(simlog.index, simCollection.index)
    push!(simlog.position, simCollection.position)
    push!(simlog.velocity, simCollection.velocity)
    push!(simlog.force, simCollection.force)

    return simCollection
end

function record_simulation_bench(simCollection, simlog)
    
    #push!(simCollection, simlog)
    #this is ultimately how the function should work, just a clean shove
    #but GenericObjectCollection may have to become a single array and not a mutable datatype
    #and then we can just push the entire array to the log
    #push!(simlog.currentstep, simCollection.currentstep)
    #@btime push!($simlog.name, $simCollection.name)
    #@btime push!($simlog.index, $simCollection.index)
    @btime append!($simlog.position, $simCollection.position)
    @btime append!($simlog.velocity, $simCollection.velocity)
    @btime append!($simlog.force, $simCollection.force)

    return simCollection
end
function record_position(positionLog, currentstep, objectname, objectindex, position)
    push!(positionLog.currentstep, currentstep)
    append!(positionLog.name, objectname)
    append!(positionLog.index, objectindex)
    append!(positionLog.position, position)
end

"""
    $unique_pairlist(a::AbstractArray)

Return a vector of static vectors, wherein each static vector is a unique pair of objects from the array of a.
In the case of position, the return pair list is a vector of static vectors of 2 mutable vectors of 3 floats each.

This is the most Naive pairlist writer.
"""
function unique_pairlist(a::AbstractArray)
    # TODO only push unique pairs to the list for eachindex(a), instead of for each pair
    list = []
    counter = 0
    j_cutoff = length(a)
    sizehint!(list, length(a)-1)

    for i in eachindex(a)
        for j in eachindex(a) 
            if j > i
                push!(list, SVector(a[i],a[j]))
            end
        end
    end
    return list
end

function generate_distance_i!(i, pairslist)
    distance_i = []
    if pairslist[1][1] == i
        for each in eachindex(pairslist)
            if pairslist[each][1] == i
                push!(distance_i, pairslist[each][3])
            elseif pairslist[each][1] !== i
                return distance_i
            end
        end
    else
        for each in eachindex(pairslist)
            if pairslist[each][1] == i
                push!(distance_i, pairslist[each][3])
            elseif pairslist[each][1] !== i
                return distance_i
            end
        end
    end
end


function force_lennardjones!(i, force,  pairslist, position)
    #TODO make epsilon and sigma user configurable 
    eps = 1
    σ = 0.0001

    #distance_i = generate_distance_i(i, pairslist)

    for each in eachindex(pairslist)
        if pairslist[each][1] == i
            d2 = pairslist[each][3]
            j = pairslist[each][2]
            # no idea if the dot is needed in .=
            d = position[i] .- position[j]

            # this is currently incorrect because d2 is the generalized distance between them, not the vectorized distances
            # we need a pairslist with the component distances not the straight-line distances
            force[i] .+= (24*eps ./ d ) .* ((2*σ ./ d).^12 .- (σ ./ d).^6)
            force[j] .-= force[i]
        end
    end
    return force
end 

function boundary_reflect!(position, velocity, collector)
    for each in eachindex(position)
    
    end  
end

function simulate!(simulation::GenericSimulation, collector)
    steps = simulation.system.duration
    
    stepwidth = simulation.system.stepwidth
    stepwidthsqrd = stepwidth^2


    logSimulation = simulation.do_logging
    simcollection = simulation.system.objectcollection

    objectcount = collector.objectnumber

    currentstep = simulation.system.objectcollection.currentstep 
    objectname = simulation.system.objectcollection.name
    objectindex = simulation.system.objectcollection.index
    mass = simulation.system.objectcollection.mass
    
    radius = simulation.system.objectcollection.radius
    position = simulation.system.objectcollection.position

    velocity = simulation.system.objectcollection.velocity
    force_currentstep = simulation.system.objectcollection.force
    force_nextstep = copy(force_currentstep)



    

    if logSimulation == true

        
        #simlog = GenericLoggerAppend(
            #currentstep,
            #objectname,
            #objectindex,
            #position,
            #velocity,
            #force_currentstep
        #)
        positionLog = PositionLogger(
            currentstep,
            objectname,
            objectindex,
            position
        )
        #pairslist = InPlaceNeighborList(x=position, cutoff=0.1, parallel=false)

        for step_n in 1:steps
            #neighborlist!(pairslist)
            pairslist = neighborlist(position, 0.1; parallel=false)
            for i in eachindex(position)

                force_lennardjones!(i, force_currentstep, pairslist, position)
                position[i] = position[i] .+ (velocity[i] .* stepwidth) .+ (force_currentstep[i] ./ mass[i] .* stepwidth^2/2)
                force_nextstep = force_currentstep .+ force_lennardjones!(i, force_currentstep, pairslist, position)
                velocity[i] = velocity[i] .+ (force_currentstep[i] .* force_nextstep[i] ./ mass[i] .* stepwidth/2)
                
            end
            #println(force_nextstep == force_currentstep) if there is at least 1 interaction, this will be false
            force_currentstep = force_nextstep
            currentstep = step_n
            record_position(positionLog, currentstep, objectname, objectindex, position)
            #@btime record_position($positionLog, $currentstep, $objectname, $objectindex, $position)
            #update!(pairslist, position)
        end
        #return positionLog
        #return simlog
    elseif logSimulation == false
        for step_n in 1:steps
            posVel_multiply!(position, velocity)
            currentstep = step_n
        end
        
    end
end


    
#end # module