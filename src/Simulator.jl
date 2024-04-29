#module Simulator
using UUIDs
using Distributions
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
mutable struct GenericLogger <: Logger 
    steparray::AbstractArray{Integer, 1}
    positionrecord::AbstractArray{AbstractArray{AbstractFloat,3}, 1}
    velocityrecord::AbstractArray{AbstractArray{AbstractFloat,3}, 1}
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

function record_simulation(step_n, position, velocity, simCollection, loopCollection, objectnumber)
    #push!(myLog,(current_step, object_index, position[1], position[2], position[3], velocity[1], velocity[2], velocity[3]))
    #c_x = position[;;; 1]
    #println(c_x)
    
    #myObjectCollection = DataFrame(
    #    current_step=step_n,
    #    object_index=collect(Int, 1:objectnumber),
    #    position_x=fill(position, length(1:objectnumber)),
    #    position_y=fill(position, length(1:objectnumber)),
    #    position_z=fill(position, length(1:objectnumber)),
    #    velocity_x=fill(velocity, length(1:objectnumber)),
    #    velocity_y=fill(velocity, length(1:objectnumber)),
    #    velocity_z=fill(velocity, length(1:objectnumber))
    #    )
    #println("Here is simCollection before cat: ", simCollection)
    #simCollection = vcat(simCollection, loopCollection)
    append!(simCollection, loopCollection)
    #println("Here is simCollection after cat: ", simCollection)
    #after cat, simcollection is right but it completely breaks when we return the value to the outerfunctions.
    #it is overridden by the previous local variable.
    #return simCollection === vcat(simCollection, loopCollection)
    return simCollection
end


function update_position!(simulation::GenericSimulation)
    position .*= velocity # this is sus. fix
    return
end

function simulate!(simulation::GenericSimulation)
    steps = simulation.system.duration
    step_n = simulation.system.currentstep
    step_size = simulation.system.stepwidth


    logSimulation = simulation.do_logging
    simCollection = simulation.system.objectcollection
    #println("simCollection is a ", typeof(simCollection))


    objectnumber = valuesCollector.objectnumber
    position = cat(simCollection.position_x, simCollection.position_y, simCollection.position_z; dims=3)
    velocity = cat(simCollection.velocity_x, simCollection.velocity_y, simCollection.velocity_z; dims=3)
    loopCollection = copy(simCollection)
    if logSimulation == true
        #myLog = simulation.system.objectcollection
        #println(myLog)
        for step_n in 1:steps
            # TODO figure out how unwrap a 3D vector and place back into a dataframe to get rid of this piecewise mess

            loopCollection.position_x .*= loopCollection.velocity_x
            loopCollection.position_y .*= loopCollection.velocity_y
            loopCollection.position_z .*= loopCollection.velocity_z
            #println("Here is the Loop collection after math: ", loopCollection)
            #position .*= velocity
            record_simulation(step_n, position, velocity, simCollection, loopCollection, objectnumber)
        end
        return simCollection
    elseif logSimulation == false
        for step_n in 1:steps
            position .*= velocity
        end
        return position, velocity
    end
end



a = [[], [], []]

    
#end # module