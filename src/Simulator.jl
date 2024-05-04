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
    currentstep::AbstractArray{Integer, 1}
    name::AbstractArray{AbstractArray{String, 1}, 1}
    #mass::AbstractArray{AbstractArray{Number, 1}, 1}
    #radius::AbstractArray{AbstractFloat, 1}
    index::AbstractArray{AbstractArray{Integer, 1}, 1}
    position::AbstractArray{AbstractArray{MVector{AbstractFloat, 3}, 1}, 1}
    velocity::AbstractArray{AbstractArray{MVector{AbstractFloat, 3}, 1}, 1}
    force::AbstractArray{bstractArray{MVector{AbstractFloat, 3}, 1}, 1}

    #uniqueID::AbstractArray{UUID,1}

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
    currentstep::Integer
    name::AbstractArray{String, 1}
    #mass::AbstractArray{Number, 1}
    #radius::AbstractArray{AbstractFloat, 1}
    index::AbstractArray{Integer, 1}
    position::AbstractArray{MVector{AbstractFloat, 3}, 1}
    velocity::AbstractArray{MVector{AbstractFloat, 3}, 1}
    force::AbstractArray{MVector{AbstractFloat, 3}, 1}

    
    #push!(simCollection, simlog)
    #this is ultimately how the function should work, just a clean shove
    #but GenericObjectCollection may have to become a single array and not a mutable datatype
    #and then we can just push the entire array to the log
    push!(simCollection.currentstep, simlog.currentstep)
    push!(simCollection.name, simlog.name)
    push!(simCollection.index, simlog.index)
    push!(simCollection.position, simlog.position)
    push!(simCollection.velocity, simlog.velocity)
    push!(simCollection.force, simlog.force)

    return simCollection
end


function update_position!(simulation::GenericSimulation)
    position .*= velocity # this is sus. fix
    return
end
function posVel_multiply!(position, velocity)
    for i in eachindex(position)
        position[i] .*= velocity[i]
    end
    return position
end

function simulate!(simulation::GenericSimulation)
    steps = simulation.system.duration
    
    step_size = simulation.system.stepwidth


    logSimulation = simulation.do_logging
    simCollection = simulation.system.objectcollection

    objectnumber = valuesCollector.objectnumber

    currentstep = simulation.system.currentstep 
    objectname = simcollection.name
    objectindex = simcollection.index
    #mass = simcollection.mass
    #radius = simCollection.radius
    position = simcollection.position
    velocity = simcollection.velocity
    force = simcollection.force

    if logSimulation == true
        simlog = GenericLogger(
            [currentstep],
            [objectname],
            [objectindex],
            [position],
            [velocity],
            [force]
        )
        for step_n in 1:steps

            posVel_multiply!(position, velocity)
            record_simulation(simCollection, simlog)
            currentstep = step_n
        end
        return simlog
    elseif logSimulation == false
        for step_n in 1:steps
            posVel_multiply!(position, velocity)

        end
        return position, velocity
    end
end


    
#end # module