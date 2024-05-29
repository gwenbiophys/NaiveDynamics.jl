#module Simulator

#println("Starting up!")
export 
    System,
    Simulation,
    TimeMethod,
    GenericIntegrator,
    Logger,
    GenericLogger,
    LoggerChunk,
    GenericSystem,
    GenericSimulation,
    simulate!,
    simulate_unified!,
    simulate_oneloop!
    #record_simulation



abstract type System end
abstract type Simulation end
abstract type TimeMethod end
struct GenericIntegrator <: TimeMethod
    stepwidth::Integer
end


abstract type Logger end

struct LoggerChunk <: Logger 
    chunk::Vector{GenericObjectCollection}
end

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
    logChunkLength::Integer
end

#function record_simulation(step_n, position, velocity, myLog)
    
 #   push!(myLog.steparray, step_n)
 #   push!(myLog.positionrecord, position)
   # push!(myLog.velocityrecord, velocity)
  #  return #println("You are logging!")
#end

function write_chunk!(simChunk)

end

function record_simulation(step_n, chunk_length, simCollection, simlog::LoggerChunk)
    chunk_index = 2
    if chunk_index != chunk_length

        chunk_index += 1
    elseif chunk_index === chunk_length
        write_chunk!(simlog)
        chunk_index = 1
    end
    
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
    eps = 0.0001
    σ = 0.0001

    #distance_i = generate_distance_i(i, pairslist)
``
    #neighborlist() fails when it has zero neighbors, this is a temporary fix
    if length(pairslist)  < 1
        return
    end

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
    #return force
end 

function boundary_reflect!(ithCoord, ithVelo, collector::Collector)
    # can this be evaluated more efficiently?
    #restructureing would allow a simple forloop
    if collector.min_xDim > ithCoord[1] 
        ithVelo[1] = -ithVelo[1] 
        ithCoord[1] = collector.min_xDim
    end
    if collector.max_xDim < ithCoord[1] 
        ithVelo[1] = -ithVelo[1] 
        ithCoord[1] = collector.max_xDim
    end

    if collector.min_yDim > ithCoord[2] 
        ithVelo[2] = -ithVelo[2] 
        ithCoord[2] = collector.min_yDim
    end
    if collector.max_yDim < ithCoord[2] 
        ithVelo[2] = -ithVelo[2] 
        ithCoord[2] = collector.max_yDim
    end

    if collector.min_zDim > ithCoord[3] 
        ithVelo[3] = -ithVelo[3] 
        ithCoord[3] = collector.min_zDim
    end
    if collector.max_zDim < ithCoord[3] 
        ithVelo[3] = -ithVelo[3] 
        ithCoord[3] = collector.max_zDim
    end
end

@inline function boundary_reflect!(position::Vec3D, velocity::Vec3D, collector::Collector)
    # can this be evaluated more efficiently?
    #restructureing would allow a simple forloop
    # this function already evaluates in 300 ns at zero allocs for 100 atoms. And it is very readily parallelizable. 
    # sooooooooooooooooooooooooo there'd have to be some substantial work to improve it
    for each in eachindex(position)
        if collector.min_xDim > position[each][1] 
            velocity[each][1] = -velocity[each][1] 
            position[each][1] = collector.min_xDim
        end
        if collector.max_xDim < position[each][1] 
            velocity[each][1] = -velocity[each][1] 
            position[each][1] = collector.max_xDim
        end

        if collector.min_yDim > position[each][2] 
            velocity[each][2] = -velocity[each][2] 
            position[each][2] = collector.min_yDim
        end
        if collector.max_yDim < position[each][2] 
            velocity[each][2] = -velocity[each][2] 
            position[each][2] = collector.max_yDim
        end

        if collector.min_zDim > position[each][3] 
            velocity[each][3] = -velocity[each][3] 
            position[each][3] = collector.min_zDim
        end
        if collector.max_zDim < position[each][3] 
            velocity[each][3] = -velocity[each][3] 
            position[each][3] = collector.max_zDim
        end
    end
end
⊗(a, b) = a .* b

function dumloop_add!(d::Vec3D, e::Vec3D)
    for i in eachindex(d)
        d[i] .+= e[i]
    end
end

function dumloop_multiply!(d::Vec3D, e::Vec3D)
    for i in eachindex(d)
        d[i] .*= e[i]
    end
end
function dumloop_multiply!(d::Vec3D, e::Vector{T}) where T
    for i in eachindex(d)
        d[i] .*= e[i]
    end
end

function dumloop_product(result::StatVec3D, d::StatVec3D, e::StatVec3D )
    for i in eachindex(d)
        result[i] = d[i] .* e[i]
    end
end
function dumloop_product(result::Vec3D, d::Vec3D, e::Vec3D )
    for i in eachindex(d)
        result[i] .= d[i] .* e[i]
    end
end
function dumloop_multiply!(d::Vec3D, e::Number)
    for i in eachindex(d)
        d[i] .*= e
    end
end


function dumloop_divide!(d, e)
    for i in eachindex(d)
        d[i] ./= e[i]
    end
  end

"""
    simulate!(simulation::GenericSimulation, collector)

Using MVectors containing dimensional data, simulate by performing calculations for each property of particles
and waiting until completion to advance to the next property/interaction.
At May 27th, best performing and no memory scaling with duration.
"""

function simulate!(simulation::GenericSimulation, collector)
    steps = simulation.system.duration
    
    stepwidth = simulation.system.stepwidth
    stepwidthsqrd = stepwidth^2



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
    inverse_mass = 1 ./ copy(mass)


    #for slice in 1:simulation.logChunkLength 

    #end

    positionIntermediate1 = copy(position)
    positionIntermediate2 = copy(position)
    velocityIntermediate1 = copy(velocity)
    for each in eachindex(positionIntermediate1)
        zero(positionIntermediate1[each])
        zero(positionIntermediate2[each])
        zero(velocityIntermediate1[each])
    end
    
    stepwidthHalf = stepwidth/2
    stepwidthSqrdHalf = stepwidth^2/2

    steps_array = zeros(Int8, steps)

    #pairslist = InPlaceNeighborList(x=position, cutoff=0.1, parallel=false)
    pairslist = []
    for step_n in eachindex(steps_array)
        #neighborlist!(pairslist)
        #try
        #pairslist = neighborlist(position, 0.02;)
        
        #catch
            #pairslist = []
        #end
        #force_lennardjones!(i, force_currentstep, pairslist, position)
        positionIntermediate1 = velocity

        dumloop_multiply!(positionIntermediate1, stepwidth)
        positionIntermediate2 = force_currentstep 
        dumloop_divide!(positionIntermediate2, mass) 
        dumloop_multiply!(positionIntermediate2, inverse_mass)
        dumloop_multiply!(positionIntermediate1, stepwidthSqrdHalf)
        dumloop_add!(position, positionIntermediate1)  
        dumloop_add!(position, positionIntermediate2)

        velocityIntermediate1[i] .= stepwidthHalf .* ((force_currentstep[i] .* force_nextstep[i]) ./ mass[i])
        velocity[i] .+= velocityIntermediate1[i]
        

        force_nextstep = force_currentstep

        #force_lennardjones!(i, force_nextstep, pairslist, position)
        dumloop_product(velocityIntermediate1, force_currentstep, force_nextstep )
        dumloop_multiply!(velocityIntermediate1, inverse_mass)
        dumloop_multiply!(velocityIntermediate1, stepwidthHalf)
        dumloop_add!(velocity, velocityIntermediate1)

        #internally, this adds zero allocations, but it adds 1 per step
        boundary_reflect!(position, velocity, collector)

        #for i in eachindex(objectindex)
        
            
            #bc force_lennardjones acts on more than 1 force object at a time, this cannot be considered threadsafe
            # without additional protections the way forward may be to allocate minicopies of JLForce
            # before entering this loop to predict their sizes. and then sum LJ forces internally before 
            #summing them to the overall forces
            #but im just writing words, no idea if that would even work anywhere
            #force_lennardjones!(i, force_currentstep, pairslist, position)

            
            #position[i] += velocity[i] * stepwidth .+ force_currentstep[i] ./ mass[i] .* stepwidth^2/2

            #force_nextstep = force_currentstep
            #force_lennardjones!(i, force_nextstep, pairslist, position)
            #velocity[i] .+= (force_currentstep[i] .* force_nextstep[i] ./ mass[i] .* stepwidth/2)

            # QUERY force should be *dumped* after each application as they are applied into the simulation?
            # the present implementation adds forces endlessly

            #boundary_reflect!(position[i], velocity[i], collector)


        #end
        #println(force_nextstep == force_currentstep) if there is at least 1 interaction, this will be false
        #force_currentstep = force_nextstep
        
        #currentstep = step_n

        #record_position(positionLog, currentstep, objectname, objectindex, position)
        #@btime record_position($positionLog, $currentstep, $objectname, $objectindex, $position)
        #update!(pairslist, position)
    end
    #return simlog

end

"""
    simulate_unified!(simulation::GenericSimulation, collector)

Using SVectors containing dimensional data, simulate by performing calculations for each property of particles
and waiting until completion to advance to the next property/interaction.
At May 27th, performance still broken.
"""

function simulate_unified!(simulation::GenericSimulation, collector)
    steps = simulation.system.duration
    
    stepwidth = simulation.system.stepwidth
    stepwidthsqrd = stepwidth^2



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

    #for slice in 1:simulation.logChunkLength 

    #end

    positionIntermediate1 = copy(position)
    positionIntermediate2 = copy(position)
    velocityIntermediate1 = copy(velocity)
    for each in eachindex(positionIntermediate1)
        zero(positionIntermediate1[each])
        zero(positionIntermediate2[each])
        zero(velocityIntermediate1[each])
    end
    
    stepwidthHalf = stepwidth/2
    stepwidthSqrdHalf = stepwidth^2/2

    steps_array = zeros(Int8, steps)

    #pairslist = InPlaceNeighborList(x=position, cutoff=0.1, parallel=false)
    pairslist = []
    for step_n in eachindex(steps_array)
        #neighborlist!(pairslist)
        #try
        #pairslist = neighborlist(position, 0.02;)
        
        #catch
            #pairslist = []
        #end
        #force_lennardjones!(i, force_currentstep, pairslist, position)
        

#copied verbatim from Molly.jl, test for allocations 
#sys.coords += sys.velocities .* sim.dt .+ ((accels_t .* sim.dt ^ 2) ./ 2)

            #accels_t = forces_t ./ masses(sys)


        #sys.velocities += ((accels_t .+ accels_t_dt) .* sim.dt / 2)


        positionIntermediate1 = velocity
        positionIntermediate1 .*= stepwidth 
        positionIntermediate2 = force_currentstep
        positionIntermediate2 ./= mass 
        positionIntermediate2 .*= stepwidthSqrdHalf 
        position .+= positionIntermediate1

        position .+= positionIntermediate2


        force_nextstep = force_currentstep
        #force_lennardjones!(i, force_nextstep, pairslist, position)
        #velocity .= (force_currentstep .⊗ force_nextstep) ./ mass * stepwidth/2

        #velocityIntermediate1 = force_currentstep .* force_nextstep
        dumloop_product(velocityIntermediate1, force_currentstep, force_nextstep)
        velocityIntermediate1 ./= mass
        velocityIntermediate1 .*= stepwidthHalf 
        velocity .+= velocityIntermediate1 

        #boundary_reflect!(position[i], velocity[i], collector)

        #for i in eachindex(objectindex)
        
            
            #bc force_lennardjones acts on more than 1 force object at a time, this cannot be considered threadsafe
            # without additional protections the way forward may be to allocate minicopies of JLForce
            # before entering this loop to predict their sizes. and then sum LJ forces internally before 
            #summing them to the overall forces
            #but im just writing words, no idea if that would even work anywhere
            #force_lennardjones!(i, force_currentstep, pairslist, position)

            
            #position[i] += velocity[i] * stepwidth .+ force_currentstep[i] ./ mass[i] .* stepwidth^2/2

            #force_nextstep = force_currentstep
            #force_lennardjones!(i, force_nextstep, pairslist, position)
            #velocity[i] .+= (force_currentstep[i] .* force_nextstep[i] ./ mass[i] .* stepwidth/2)

            # QUERY force should be *dumped* after each application as they are applied into the simulation?
            # the present implementation adds forces endlessly

            #boundary_reflect!(position[i], velocity[i], collector)


        #end
        #println(force_nextstep == force_currentstep) if there is at least 1 interaction, this will be false
        #force_currentstep = force_nextstep
        
        #currentstep = step_n

        #record_position(positionLog, currentstep, objectname, objectindex, position)
        #@btime record_position($positionLog, $currentstep, $objectname, $objectindex, $position)
        #update!(pairslist, position)
    end
    #return simlog

end

"""
    simulate_oneloop!(simulation::GenericSimulation, collector)

Using MVectors as the containers of dimensional data, simulate by assigning a block of calculations to each particle. 
At May 27th, extremely allocation heavy and slow.

"""

function simulate_oneloop!(simulation::GenericSimulation, collector)
    steps = simulation.system.duration
    
    stepwidth = simulation.system.stepwidth
    stepwidthsqrd = stepwidth^2



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

    #for slice in 1:simulation.logChunkLength 

    #end

    positionIntermediate1 = copy(position)
    positionIntermediate2 = copy(position)
    velocityIntermediate1 = copy(velocity)
    for each in eachindex(positionIntermediate1)
        zero(positionIntermediate1[each])
        zero(positionIntermediate2[each])
        zero(velocityIntermediate1[each])
    end
    
    stepwidthHalf = stepwidth/2
    stepwidthSqrdHalf = stepwidth^2/2

    steps_array = zeros(Int8, steps)

    #pairslist = InPlaceNeighborList(x=position, cutoff=0.1, parallel=false)
    pairslist = []
    for step_n in eachindex(steps_array)
        #neighborlist!(pairslist)
        #try
        #pairslist = neighborlist(position, 0.02;)
        
        #catch
            #pairslist = []
        #end


        #boundary_reflect!(position[i], velocity[i], collector)

        for i in eachindex(objectindex)
        
            
            #bc force_lennardjones acts on more than 1 force object at a time, this cannot be considered threadsafe
            # without additional protections the way forward may be to allocate minicopies of JLForce
            # before entering this loop to predict their sizes. and then sum LJ forces internally before 
            #summing them to the overall forces
            #but im just writing words, no idea if that would even work anywhere
            #force_lennardjones!(i, force_currentstep, pairslist, position)
            

            #positionIntermediate1[i] = velocity[i] 
            #positionIntermediate1[i] .*= stepwidth
            #positionIntermediate2[i] = force_currentstep[i] 
            #positionIntermediate2[i] ./= mass[i] 
            #positionIntermediate2[i] .*= stepwidthSqrdHalf
            #position[i] .+= positionIntermediate1[i] 
            #position[i] .+= positionIntermediate2[i]
            positionIntermediate1[i] .= stepwidth .* velocity[i]
            positionIntermediate2[i] .= stepwidthSqrdHalf .* (force_currentstep[i] ./ mass[i])
            positionIntermediate1[i] .+= positionIntermediate2[i]
            position[i] .+= positionIntermediate1[i]

            

            #force_nextstep[i] = force_currentstep[i]
            #force_lennardjones!(i, force_nextstep, pairslist, position)
            
            #velocityIntermediate1[i] .= force_currentstep[i] .* force_nextstep[i]
            #velocityIntermediate1[i] ./= mass[i]
            #velocityIntermediate1[i] .*= stepwidthHalf
            #velocity[i] .+= velocityIntermediate1[i]
            velocityIntermediate1[i] .= stepwidthHalf .* ((force_currentstep[i] .* force_nextstep[i]) ./ mass[i])
            velocity[i] .+= velocityIntermediate1[i]

            # QUERY force should be *dumped* after each application as they are applied into the simulation?
            # the present implementation adds forces endlessly

            #boundary_reflect!(position[i], velocity[i], collector)


        end
        #println(force_nextstep == force_currentstep) if there is at least 1 interaction, this will be false
        #force_currentstep = force_nextstep
        
        #currentstep = step_n

        #record_position(positionLog, currentstep, objectname, objectindex, position)
        #@btime record_position($positionLog, $currentstep, $objectname, $objectindex, $position)
        #update!(pairslist, position)
    end
    #return simlog

end
    
#end # module