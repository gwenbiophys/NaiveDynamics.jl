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

push!
function update_position!(simulation::GenericSimulation)
    position .*= velocity # this is sus. fix
    return
end
function posVel_multiply!(position, velocity)
    #println(position)
    #println(typeof(position))
    #println()
    #position = position .* velocity'
    #for i in 1:length(position)
        #position[i] *= velocity[i]
    #end
    return position
end

function simulate!(simulation::GenericSimulation, collector)
    steps = simulation.system.duration
    
    stepwidth = simulation.system.stepwidth
    stepwidthsqrd = stepwidth^2


    logSimulation = simulation.do_logging
    simcollection = simulation.system.objectcollection

    objectnumber = collector.objectnumber

    currentstep = simulation.system.objectcollection.currentstep 
    objectname = simulation.system.objectcollection.name
    objectindex = simulation.system.objectcollection.index
    mass = simulation.system.objectcollection.mass
    
    #radius = simCollection.radius
    position = simulation.system.objectcollection.position

    velocity = simulation.system.objectcollection.velocity
    force_currentstep = simulation.system.objectcollection.force

    if logSimulation == true

        #simlog = GenericLogger(
        #    [currentstep],
         #   [objectname],
         #   [objectindex],
         #   [position],
         #   [velocity],
         #   [force]
        #)
        simlog = GenericLoggerAppend(
            currentstep,
            objectname,
            objectindex,
            position,
            velocity,
            force_currentstep
        )

        for step_n in 1:steps
            for i in eachindex(position)
                #force_currentstep calculations and sum of force calculations
                force_position = force_currentstep[i] ./ mass[i] .* stepwidth^2/2 #can this be combined into next line, and why would i do that TO REDUCE COPYING
                position[i] = position[i] .+ (velocity[i] .* stepwidth) .+ (force_currentstep[i] ./ mass[i] .* stepwidth^2/2)
                #force_nextstep= sum of new applied energies
                force_nextstep = force_currentstep
                velocity[i] = velocity[i] .+ (force_currentstep[i] .* force_nextstep[i] ./ mass[i] .* stepwidth/2)
                currentstep = step_n
            end
        end
        #return simlog
    elseif logSimulation == false
        for step_n in 1:steps
            posVel_multiply!(position, velocity)
            currentstep = step_n
        end
        
    end
end


    
#end # module