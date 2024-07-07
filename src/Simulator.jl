#module Simulator

#println("Starting up!")
export 

    SimulationSpecification,
    GenericSpec,
    simulate!,
    simulate_bravado!,
    simulate_fusedBravado!,
    simulate_MVec!,
    simulate_SVec!,
    simulate_oneloop!#,
    #record_simulation,
    #update_chunk!,
    #write_chunk!



abstract type SimulationSpecification end
struct GenericSpec <: SimulationSpecification
    duration::Integer
    stepwidth::Integer
    currentstep::Integer
    logChunkLength::Integer
end


function write_chunk!(simLog, simChunk)
    return append!(simLog, simChunk)
end
function update_chunk!(chunk_index, simChunk, objectCollection)
    #return simChunk[chunk_index] = objectCollection
    return fill!(simChunk[chunk_index], objectCollection)
end

function record_simulation(step_n, chunk_index, chunk_length, simChunk, simLog, objectCollection)

    if chunk_length > chunk_index
        simChunk = update_chunk!(chunk_index, simChunk, objectCollection)
        #println(simChunk[5].position)
        return chunk_index += 1 
    elseif chunk_index == chunk_length
        simLog = write_chunk!(simLog, simChunk)
        #println(simLog[each].position) for each in 1:length(simLog)
        return chunk_index = 1
    else
        println("something got bungle grundled in record_simulation")
    end


end

function record_position(positionLog, currentstep, objectname, objectindex, position)
end

"""
    unique_pairlist(a::AbstractArray)

Return a vector of tuples, wherein each tuple 
contains the indices and squared distance of pairs within a threshold float and the component distances.
In the case of position, the return pair list is a vector of static vectors of 2 mutable vectors of 3 floats each.

This is the most Naive pairlist writer.
"""
function unique_pairlist!(a::AbstractArray, threshold::AbstractFloat)
    # TODO only push unique pairs to the list for eachindex(a), instead of for each pair
    list = []
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
                    push!(list, tuple(i, j, dx, dy, dz, d2)) #could have pairlist be arbitarily large and just set to zero? eh. patience.
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


function force_lennardjones!(force::Vec3D,  pairslist, position)
    #TODO make epsilon and sigma user configurable 
    eps = 1.0f-23
    #eps = 0.00000001 any value diff from this causes instant chaos
    σ = 0.3

    #distance_i = generate_distance_i(i, pairslist)

    #neighborlist() fails when it has zero neighbors, this is a temporary fix
    if length(pairslist) < 1
        return
    end

    for each in eachindex(pairslist)
    
        #d2 = pairslist[each][3]
        i = pairslist[each][1]
        j = pairslist[each][2]

        d = position[i] .- position[j]


        force[i] .+= (24*eps ./ d ) .* ((2*σ ./ d).^12 .- (σ ./ d).^6)
        force[j] .-= force[i]
        
    end
    #println(typeof(force))
    #return force
end

function sum_forces!(force, force1)
    dumloop_add!(force, force1)

    for each in eachindex(force1)
        zero(force1[each])
    end
    return force
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

function boundary_reflect!(position::Vec3D, velocity::Vec3D, collector::GenericRandomCollector)
    # does not actually reflect, just reverses the particle
    for each in eachindex(position)
        if collector.min_xDim > position[each][1] 
            velocity[each][1] *= -1
            position[each][1] = collector.min_xDim
        end
        if collector.max_xDim < position[each][1] 
            velocity[each][1] *= -1
            position[each][1] = collector.max_xDim
        end

        if collector.min_yDim > position[each][2] 
            velocity[each][2] *= -1 
            position[each][2] = collector.min_yDim
        end
        if collector.max_yDim < position[each][2] 
            velocity[each][2] *= -1
            position[each][2] = collector.max_yDim
        end

        if collector.min_zDim > position[each][3] 
            velocity[each][3] *= -1 
            position[each][3] = collector.min_zDim
        end
        if collector.max_zDim < position[each][3] 
            velocity[each][3] *= -1
            position[each][3] = collector.max_zDim
        end
    end
end

"""
    simulate!(simulation::GenericSimulation, collector)

Using MVectors containing dimensional data, simulate by performing calculations for each property of particles
and waiting until completion to advance to the next property/interaction.
At June 29th, best performing and minimal allocations per time step. better performing for SVecs 
but does not interoperate with other functions like boundary_reflect!() or forces.
"""

function simulate!(simulation, collector)


    steps::Int64 = simulation.system.duration
    
    stepwidth::Int64 = simulation.system.stepwidth
    

    #objectcount = collector.objectnumber

    #currentstep = simulation.system.objectcollection.currentstep 
    #objectname = simulation.system.objectcollection.name
    #objectindex = simulation.system.objectcollection.index
    mass::Vector{Int64} = simulation.system.objectcollection.mass
    
    #radius = simulation.system.objectcollection.radius
    stepwidthHalf::Float32 = stepwidth/2
    stepwidthSqrdHalf::Float32 = stepwidth^2/2

    position::Vec3D{Float32} = simulation.system.objectcollection.position

    velocity::Vec3D{Float32} = simulation.system.objectcollection.velocity
    #velocity = @views simulation.system.objectcollection.velocity
    force_currentstep::Vec3D{Float32} = simulation.system.objectcollection.force
    force_nextstep::Vec3D{Float32} = copy(force_currentstep)
    force_LJ::Vec3D{Float32} = copy(force_currentstep)
    inverse_mass::Vector{Float32} = 1 ./ mass


    chunk_index::Int64 = 2
    chunk_length::Int64 = 10
    simLog::Vector{GenericObjectCollection} = []
    simChunk::Vector{GenericObjectCollection} = [simulation.system.objectcollection for _ in 1:chunk_length]
    
    
    #for slice in 1:simulation.logChunkLength 

    #end

    positionIntermediate1::Vec3D{Float32} = copy(position)
    positionIntermediate2::Vec3D{Float32} = copy(position)
    velocityIntermediate1::Vec3D{Float32} = copy(velocity)
    for each in eachindex(positionIntermediate1)
        zero(positionIntermediate1[each])
        zero(positionIntermediate2[each])
        zero(velocityIntermediate1[each])
    end

    
    #pairslist = InPlaceNeighborList(x=position, cutoff=0.1, parallel=false)
    pairslist = []
    for step_n in 1:steps

        #neighborlist!(pairslist)

        #pairslist = neighborlist(position, 0.02;)
        pairslist = unique_pairlist!(position, 0.05)

        force_lennardjones!(force_LJ, pairslist, position)
        sum_forces!(force_currentstep, force_LJ)

        positionIntermediate1 = copy(velocity) # TODO fix

        dumloop_multiply!(positionIntermediate1, stepwidth)
        positionIntermediate2 = force_currentstep 
        dumloop_divide!(positionIntermediate2, mass) 
        dumloop_multiply!(positionIntermediate2, inverse_mass)
        dumloop_multiply!(positionIntermediate1, stepwidthSqrdHalf)
        dumloop_add!(position, positionIntermediate1)  
        dumloop_add!(position, positionIntermediate2)

        force_nextstep = force_currentstep


        force_lennardjones!(force_LJ, pairslist, position)
        sum_forces!(force_nextstep, force_LJ)


        dumloop_product(velocityIntermediate1, force_currentstep, force_nextstep )
        dumloop_multiply!(velocityIntermediate1, inverse_mass)
        dumloop_multiply!(velocityIntermediate1, stepwidthHalf)
        dumloop_add!(velocity, velocityIntermediate1)

        
        boundary_reflect!(position, velocity, collector)



       
        #println(force_nextstep == force_currentstep) if there is at least 1 interaction, this will be false
        force_currentstep = force_nextstep
        
        currentstep = step_n
        

        chunk_index = record_simulation(step_n, chunk_index, chunk_length, simChunk, simLog, simulation.system.objectcollection)

        #record_position(positionLog, currentstep, objectname, objectindex, position)
        #@btime record_position($positionLog, $currentstep, $objectname, $objectindex, $position)
        #update!(pairslist, position)
    end
    return simLog
end

function simulate_bravado!(sys::GenericObjectCollection, spec::GenericSpec, clct::GenericRandomCollector)

    stepwidthHalf::Float32 = spec.stepwidth/2
    stepwidthSqrdHalf::Float32 = spec.stepwidth^2/2

    force_nextstep::Vec3D{Float32} = copy(sys.force)
    force_LJ::Vec3D{Float32} = copy(sys.force)
    inverse_mass::Vector{Float32} = 1 ./ sys.mass


    chunk_index::Int64 = 2

    simLog::Vector{GenericObjectCollection} = []
    simChunk::Vector{GenericObjectCollection} = [sys for _ in 1:spec.logChunkLength]
    
    poslog::Vector{Vec3D{Float32}} = []
    push!(poslog, copy.(sys.position))
    positionIntermediate1::Vec3D{Float32} = copy.(sys.velocity)
    positionIntermediate2 = copy.(sys.position)

    velocityIntermediate1::Vec3D{Float32} = copy.(sys.velocity)

    for each in eachindex(positionIntermediate2)
        zero(positionIntermediate2[each])
        zero(velocityIntermediate1[each])
    end

    
    
    #pairslist = InPlaceNeighborList(x=position, cutoff=0.1, parallel=false)
    pairslist = []

    
    for step_n in 1:spec.duration

        #neighborlist!(pairslist)

        #pairslist = neighborlist(position, 0.02;)
        pairslist = unique_pairlist!(sys.position, 0.3)

        force_lennardjones!(force_LJ, pairslist, sys.position)
        sys.force = sum_forces!(sys.force, force_LJ)

        positionIntermediate1 = copy.(sys.velocity)
        #fill!(positionIntermediate1, sys.velocity)
        #println(sys.position)
        #println(sys.velocity)
        dumloop_multiply!(positionIntermediate1, spec.stepwidth)
        positionIntermediate2 = copy.(sys.force) 
        #dumloop_divide!(positionIntermediate2, sys.mass) 
        dumloop_multiply!(positionIntermediate2, inverse_mass)
        dumloop_multiply!(positionIntermediate1, stepwidthSqrdHalf)
        #println(positionIntermediate1)
        #println(sys.position)
        #println()
        dumloop_add!(sys.position, positionIntermediate1)
       # println(positionIntermediate1)
        #println(sys.position)  
        dumloop_add!(sys.position, positionIntermediate2)
        #println("velocities at simulate!() for each step")
        #println(sys.velocity)
        #println()
        force_nextstep = copy.(sys.force)


        force_lennardjones!(force_LJ, pairslist, sys.position)
        force_nextstep = sum_forces!(force_nextstep, force_LJ)

        #println("velocity before we intermediate???")
        #println(sys.velocity)
        dumloop_product(velocityIntermediate1, sys.force, force_nextstep )
        #println("velocity before we intermediate???")
        #println(sys.velocity)
        dumloop_multiply!(velocityIntermediate1, inverse_mass)
        dumloop_multiply!(velocityIntermediate1, stepwidthHalf)

        dumloop_add!(sys.velocity, velocityIntermediate1)

        
        boundary_reflect!(sys.position, sys.velocity, clct)



        #println(sys.velocity[1] == sys.velocity[2] )
        #println(force_nextstep == sys.force) #if there is at least 1 interaction, this will be false
       
            #sys.force = force_nextstep
        #for each in eachindex(sys.force)
           # sys.force[each] = force_nextstep[each]
        #end
        sys.force = copy.(force_nextstep)

        # TODO this needs to be ironed out
        push!(sys.currentstep, step_n)
        

        #chunk_index = record_simulation(step_n, chunk_index, spec.logChunkLength, simChunk, simLog, sys)
        
        #push!(simLog, sys)
        
        push!(poslog, copy.(sys.position))


        #record_position(positionLog, currentstep, objectname, objectindex, position)
        #@btime record_position($positionLog, $currentstep, $objectname, $objectindex, $position)
        #update!(pairslist, position)
    end
    println(sys.position)


    #return simLog
    return poslog
end
function simulate_fusedBravado!(sys::GenericObjectCollection, spec::GenericSpec, clct::GenericRandomCollector)

    stepwidthHalf::Float32 = spec.stepwidth/2
    stepwidthSqrdHalf::Float32 = spec.stepwidth^2/2

    force_nextstep::Vec3D{Float32} = copy(sys.force)
    force_LJ::Vec3D{Float32} = copy(sys.force)
    inverse_mass::Vector{Float32} = 1 ./ sys.mass


    chunk_index::Int64 = 2

    simLog::Vector{GenericObjectCollection} = []
    simChunk::Vector{GenericObjectCollection} = [sys for _ in 1:spec.logChunkLength]
    

    positionIntermediate1::Vec3D{Float32} = copy(sys.position)
    positionIntermediate2::Vec3D{Float32} = copy(sys.position)
    velocityIntermediate1::Vec3D{Float32} = copy(sys.velocity)
    accels_t::Vec3D{Float32} = copy(sys.force)
    accels_t_dt::Vec3D{Float32} = copy(sys.force)
    for each in eachindex(positionIntermediate1)
        zero(positionIntermediate1[each])
        zero(positionIntermediate2[each])
        zero(velocityIntermediate1[each])
        zero(accels_t[each])
        zero(accels_t_dt[each])
    end

    
    
    #pairslist = InPlaceNeighborList(x=position, cutoff=0.1, parallel=false)
    pairslist = []
    for step_n in 1:spec.duration

        #neighborlist!(pairslist)

        #pairslist = neighborlist(position, 0.02;)
        pairslist = unique_pairlist!(sys.position, 0.3)

        force_lennardjones!(force_LJ, pairslist, sys.position)
        sum_forces!(sys.force, force_LJ)


        (accels_t[i] = sys.force[i] ./ sys.mass[i] for i in eachindex(sys.force))

        (sys.position[i] += sys.velocity[i] .* spec.stepwidth[i] .+ ((accels_t[i] .* spec.stepwidth ^ 2) ./ 2) for i in eachindex(sys.position))

        (accels_t_dt[i] = force_nextstep[i] ./ sys.mass[i] for i in eachindex(force_nextstep))
        
        
        force_nextstep = sys.force


        force_lennardjones!(force_LJ, pairslist, sys.position)
        sum_forces!(force_nextstep, force_LJ)

        (sys.velocity[i] += ((accels_t[i] .+ accels_t_dt)[i] .* stepwidth / 2) for i in eachindex(sys.velocity))

        boundary_reflect!(sys.position, sys.velocity, clct)



       
       
            #sys.force = force_nextstep
        for each in eachindex(sys.force)
            sys.force[each] = force_nextstep[each]
        end
        
        currentstep = step_n
        

        chunk_index = record_simulation(step_n, chunk_index, spec.logChunkLength, simChunk, simLog, sys)

        #record_position(positionLog, currentstep, objectname, objectindex, position)
        #@btime record_position($positionLog, $currentstep, $objectname, $objectindex, $position)
        #update!(pairslist, position)
    end
    return simLog
end

function simulate_MVec!(simulation, collector)


    steps::Int64 = simulation.system.duration
    
    stepwidth::Int64 = simulation.system.stepwidth
    

    #objectcount = collector.objectnumber

    #currentstep = simulation.system.objectcollection.currentstep 
    #objectname = simulation.system.objectcollection.name
    #objectindex = simulation.system.objectcollection.index
    mass::Vector{Int64} = simulation.system.objectcollection.mass
    
    #radius = simulation.system.objectcollection.radius
    stepwidthHalf::Float32 = stepwidth/2
    stepwidthSqrdHalf::Float32 = stepwidth^2/2

    position::Vec3D{Float32} = simulation.system.objectcollection.position

    velocity::Vec3D{Float32} = simulation.system.objectcollection.velocity
    #velocity = @views simulation.system.objectcollection.velocity
    force_currentstep::Vec3D{Float32} = simulation.system.objectcollection.force
    force_nextstep::Vec3D{Float32} = copy(force_currentstep)
    force_LJ::Vec3D{Float32} = copy(force_currentstep)
    inverse_mass::Vector{Float32} = 1 ./ mass

    #simLog = DataFrame(currentstep, positionX, positionY, positionZ)
    chunk_length = 10
    simLog::Vector{GenericObjectCollection} = []
    simChunk::Vector{GenericObjectCollection} = [simulation.system.objectcollection for _ in 1:chunk_length]
    
    
    #for slice in 1:simulation.logChunkLength 

    #end

    positionIntermediate1::Vec3D{Float32} = copy(position)
    positionIntermediate2::Vec3D{Float32} = copy(position)
    velocityIntermediate1::Vec3D{Float32} = copy(velocity)
    for each in eachindex(positionIntermediate1)
        zero(positionIntermediate1[each])
        zero(positionIntermediate2[each])
        zero(velocityIntermediate1[each])
    end

    
    
    #pairslist = InPlaceNeighborList(x=position, cutoff=0.1, parallel=false)
    pairslist = []
    for step_n in eachindex(steps_array)
        #neighborlist!(pairslist)
        #try
        #pairslist = neighborlist(position, 0.02;)
        pairslist = unique_pairlist(position, 0.02)
        
        #catch
            #pairslist = []
        #end
        #force_lennardjones!(i, force_currentstep, pairslist, position)

        (accels_t[i] = force_currentstep[i] ./ mass[i] for i in eachindex(force_currentstep))

        
        (position[i] += velocity[i] .* stepwidth[i] .+ ((accels_t[i] .* stepwidth ^ 2) ./ 2) for i in eachindex(position))

        (accels_t_dt[i] = force_nextstep[i] ./ mass[i] for i in eachindex(force_nextstep))
        #force_nextstep = force_currentstep
        #force_lennardjones!(i, force_nextstep, pairslist, position)


        (velocity[i] += ((accels_t[i] .+ accels_t_dt)[i] .* stepwidth / 2) for i in eachindex(velocity))


        boundary_reflect!(position, velocity, collector)

       
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
    simulate_SVec!(simulation::GenericSimulation, collector)

Using SVectors containing dimensional data, simulate by performing calculations for each property of particles
and waiting until completion to advance to the next property/interaction.
At May 27th, performance still broken.
"""

function simulate_SVec!(simulation, collector)
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

    accels_t::StatVec3D = copy(force_currentstep)
    accels_t_dt::StatVec3D = copy(force_nextstep)
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

        #accels_update!(accels_t, force_currentstep, mass )
        (accels_t[i] = force_currentstep[i] ./ mass[i] for i in eachindex(force_currentstep))

        #position_update!(position, velocity, accels_t, stepwidth)
        (position[i] += velocity[i] .* stepwidth[i] .+ ((accels_t[i] .* stepwidth ^ 2) ./ 2) for i in eachindex(position))

        (accels_t_dt[i] = force_nextstep[i] ./ mass[i] for i in eachindex(force_nextstep))
        #force_nextstep = force_currentstep
        #force_lennardjones!(i, force_nextstep, pairslist, position)


        (velocity[i] += ((accels_t[i] .+ accels_t_dt)[i] .* stepwidth / 2) for i in eachindex(velocity))
 

        #boundary_reflect!(position[i], velocity[i], collector)

        
        #force_currentstep = force_nextstep
        
        #currentstep = step_n

        #record_position(positionLog, currentstep, objectname, objectindex, position)

        #update!(pairslist, position)
    end
    #return simlog

end

"""
    simulate_oneloop!(simulation::GenericSimulation, collector)

Using MVectors as the containers of dimensional data, simulate by assigning a block of calculations to each particle. 
At May 27th, extremely allocation heavy and slow.

"""

function simulate_oneloop!(simulation, collector)
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