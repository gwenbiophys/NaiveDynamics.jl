#module Simulator
using Revise
using UUIDs
using Distributions
println("Starting up!")
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
    objectcollection::GenericObjectCollection
end
struct GenericSimulation <: Simulation
    system::GenericSystem
    do_logging::Bool
end

function record_simulation(simulation::GenericSimulation)
    
    push!(myLog.steparray, step_n)
    push!(myLog.positionrecord, position)
    push!(myLog.velocityrecord, velocity)
    return println("You are logging!")
end

function update_position!(simulation::GenericSimulation)
    position .*= velocity # this is sus. fix
    return println("You are updating positions!")
end

function simulate!(simulation::GenericSimulation)
    steps = simulation.system.duration
    step_n = simulation.system.currentstep
    step_size = simulation.system.stepwidth

    position = simulation.system.objectcollection.position
    velocity = simulation.system.objectcollection.position

    logSimulation = simulation.do_logging

    if logSimulation == true
        myLog = GenericLogger([step_n], [position], [velocity])
        for step_n in steps
            
            position .*= velocity
            #record_simulation(simulation) # there is a scoping problem here, 
            #this function is not accepting scope from simulate! in an expected way
            push!(myLog.steparray, step_n)
            push!(myLog.positionrecord, position)
            push!(myLog.velocityrecord, velocity)
            #rowan says, throw the local variables to the function call, rahter than just being maxspeed
            # that ir doesnt work
        end
    elseif logSimulation == false
        position .*= velocity
    end
    return println(steps, step_n, step_size, position, velocity)
end





    
#end # module