#module Simulator

#println("Starting up!")
export 

    SimulationSpecification,
    GenericSpec,
    simulate!,
    simulate_dumloop!,
    simulate_SVec!
    #record_simulation,
    #update_chunk!,
    #write_chunk!



abstract type SimulationSpecification end
struct GenericSpec{T, K} <: SimulationSpecification
    duration::T
    stepwidth::T
    currentstep::T
    logChunkLength::T
    velocityDampening::K #idk a better way to handle this
end
function GenericSpec(;
                    inttype=Int64,
                    floattype=Float32,
                    duration,
                    stepwidth,
                    currentstep,
                    logLength=10,
                    vDamp
                        )
    return GenericSpec{inttype, floattype}(duration, stepwidth, currentstep, logLength, vDamp)
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

function update_pairslist!(a::Vec3D{T}, list) where T

    
    if length(a) == 2

        dx = a[1][1] - a[2][1]
        dy = a[1][2] - a[2][2] 
        dz = a[1][3] - a[2][3]
        d2 = sqrt(dx^2 + dy^2 + dz^2) 

        list[1] = tuple(1, 2, dx, dy, dz, d2)

    else


        for each in eachindex(list)

            i = list[each][1]
            j = list[each][2]

            xi = a[i][1]
            xj = a[j][1]

            yi = a[i][2]
            yj = a[j][2]

            zi = a[i][3]
            zj = a[j][3]

            dx = xi - xj
            dy = yi - yj
            dz = zi - zj
            d2 = sqrt(dx^2 + dy^2 + dz^2)  
            result = tuple(i, j, dx, dy, dz, d2)

            list[i] = result
           


        end
    end

    return list
end

function unique_pairs(a::Vec3D{T}) where T
    list_length = convert(Int64, (length(a)-1) * length(a) / 2)


    list = [tuple(i, j, a[1][1], a[1][1], a[1][1], a[1][1]) for i in 1:length(a)-1 for j in i+1:length(a)]
    update_pairslist!(a, list)

    return list
end


function threshold_pairs(list, threshold::T) where T
    
    return [list[i] for i in eachindex(list) if list[i][6] ≤ threshold]

end

function threshold_pairs_old(list, threshold::T) where T
    thresh_list::Vector{Tuple{Int64, Int64, T, T, T, T}} = [] # this is a silly fix
    # would it more perf-efficient to define a threshold list as long as the unique pairs list
    # at small n particles, and just reorder the threshlist between valid and invalid values
    # and jsut instruct functions to use the 'valid' region of the array?
    thresh_list = []
    for i in eachindex(list)
        # replace with named tuple?
        if list[i][6] ≤ threshold
            push!(thresh_list, list[i])
        end
    end

    return thresh_list

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


function force_lennardjones!(force::Vec3D{T},  pairslist, position) where T
    #TODO make epsilon and sigma user configurable 
    eps = -1f10

    σ = 0.0001   
    

    # this is silly
    for each in eachindex(force)
        map!(x->x, force[each], MVector{3, T}(0.0, 0.0, 0.0))
    end
    #zero(force)
    #for each in eachindex(force)
        #zero.(force[each])
    #end

    

    #neighborlist() fails when it has zero neighbors, this is a temporary fix
    if length(pairslist) < 1
        return
    end

    for each in eachindex(pairslist)
        #d2 = pairslist[each][3]
        i = pairslist[each][1]
        j = pairslist[each][2]
        dx = pairslist[each][3]
        dy = pairslist[each][4]
        dz = pairslist[each][5]

        d = pairslist[each][6]

        # i hate this, profoundly
        force[i][1] += (24*eps / dx ) * ((2*σ / dx)^12 - (σ / dx)^6)
        force[i][2] += (24*eps / dy ) * ((2*σ / dy)^12 - (σ / dy)^6)
        force[i][3] += (24*eps / dz ) * ((2*σ / dz)^12 - (σ / dz)^6)

        force[j][1] -= force[i][1]
        force[j][2] -= force[i][2]
        force[j][3] -= force[i][3]
        



        # incorrect, overall force is being applied to each component
        #force[i] .+= (24*eps ./ d ) .* ((2*σ ./ d).^12 .- (σ ./ d).^6)
        #force[j] .-= force[i]
        
    end


    #return force
end

function force_coulomb!(force::Vec3D{T}, pairslist, charge) where T
    k = 1 # coulomb constant
    for each in eachindex(force)
        map!(x->x, force[each], MVector{3, T}(0.0, 0.0, 0.0))
    end

    for each in eachindex(pairslist)
        i = pairslist[each][1]
        j = pairslist[each][2]
        dx = pairslist[each][3]
        dy = pairslist[each][4]
        dz = pairslist[each][5]

        #d = pairslist[each][6]

        force[i][1] += k * charge[i] * charge[j] / dx^2
        force[i][2] += k * charge[i] * charge[j] / dy^2
        force[i][3] += k * charge[i] * charge[j] / dz^2
        force[j][1] -= force[i][1]
        force[j][2] -= force[i][2]
        force[j][3] -= force[i][3]

    end
end

function sum_forces!(force, force1, force2)

    for each in eachindex(force)
        force[each] .= force1[each] .+ force2[each] 
    end

    return force
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
    rescale_velocity(velocity, Tf, gamma)

For gamma=1, we will achieve full rescaling of velocity at the end of each step.
No rescaling for zero.
"""

function rescale_velocity!(velocity::Vec3D{T}, Tf::T, γ::T, mass::Vector{T}, objectcount::Int64) where T
    Ti::T = 0.0
    kb::T = 1.0
    v::T = 0.0

    objects = convert(T, objectcount)

    
    for each in eachindex(velocity)

        v = (velocity[each][1]^2 + velocity[each][2]^2 + velocity[each][3]^2) ^ 0.5
        #println("here is v ", v)

        Ti += (2/(3 * objects * kb)) * v * mass[each]/2
        #println("Ti after calc ", Ti)
    end

    β = (1 + γ * (Tf/Ti - 1) ) ^ 0.5
    #β = (Tf/Ti) ^ 0.5
    #println("here is β ", β)

    for each in eachindex(velocity)
        velocity[each] .*= β
    end

end

"""
    simulate!(simulation::GenericSimulation, collector)

Using MVectors containing dimensional data, simulate by performing calculations for each property of particles
and waiting until completion to advance to the next property/interaction.
At June 29th, best performing and minimal allocations per time step. better performing for SVecs 
but does not interoperate with other functions like boundary_reflect!() or forces.
"""

function simulate!(sys::GenericObjectCollection, spec::GenericSpec, clct::GenericRandomCollector{T}) where T


    force_nextstep = deepcopy(sys.force)::Vec3D{T}
    force_LJ = deepcopy(sys.force)::Vec3D{T}
    force_C = deepcopy(sys.force)::Vec3D{T}

    chunk_index::Int64 = 2

# type assert error, have to decide if simLog is instantiated with a sys or not
    #simLog = []::Vector{GenericObjectCollection}
   # simChunk = [sys for _ in 1:spec.logChunkLength]::Vector{GenericObjectCollection}
    
    poslog = [sys.position]::Vector{Vec3D{T}}
    #poslog = [sys.position for each in 1:spec.duration]
    #sizehint!(poslog, spec.duration)
    #push!(poslog, copy.(sys.position)::Vector{Vec3D{Float32}})
    accels_t = copy.(sys.force)::Vec3D{T}
    accels_t_dt = copy.(sys.force)::Vec3D{T}

    #pairslist = InPlaceNeighborList(x=position, cutoff=0.1, parallel=false)

    pairslist = unique_pairs(sys.position)


    for step_n in 1:spec.duration
        

        #neighborlist!(pairslist)

        #pairslist = neighborlist(position, 0.02;)

        LJ_pairs = threshold_pairs(pairslist, convert(T, 1.0))


        force_lennardjones!(force_LJ, LJ_pairs, sys.position)
        force_coulomb!(force_C, pairslist, sys.charge)

        sum_forces!(sys.force, force_LJ, force_C)
        

        for i in eachindex(accels_t)
            accels_t[i] .= sys.force[i] ./ sys.mass[i]
        end

        for i in eachindex(sys.position)
            sys.position[i] .+= sys.velocity[i] .* spec.stepwidth .+ ((accels_t[i] .* spec.stepwidth ^ 2) ./ 2)
        end

        for i in eachindex(sys.position)
            accels_t_dt[i] .= force_nextstep[i] ./ sys.mass[i]
        end
        
        map!(x->x, force_nextstep, sys.force)



        force_lennardjones!(force_LJ, pairslist, sys.position)
        force_coulomb!(force_C, pairslist, sys.charge)

        sum_forces!(force_nextstep, force_LJ, force_C)


        for i in eachindex(sys.velocity)
            sys.velocity[i] .+= (accels_t[i] .+ accels_t_dt[i]) .* spec.stepwidth / 2
        end
        boundary_reflect!(sys.position, sys.velocity, clct)

        map!(x->x, sys.force, force_nextstep)



        currentstep = 1:step_n



        #### end of step processing

        #update!(pairslist, position)



        update_pairslist!(sys.position, pairslist)

        if step_n % 10 == 0
            rescale_velocity!(sys.velocity, clct.temperature, spec.velocityDampening, sys.mass, clct.objectnumber)
        end

        push!(poslog, deepcopy(sys.position))
        #copyto!(poslog[step_n], sys.position) # this implementation is broken atm
        #chunk_index = record_simulation(step_n, chunk_index, spec.logChunkLength, simChunk, simLog, sys)
    end

    return poslog
end


function simulate_dumloop!(sys::GenericObjectCollection, spec::GenericSpec, clct::GenericRandomCollector)

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
        #sys.force = sum_forces!(sys.force, force_LJ)

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
        #force_nextstep = sum_forces!(force_nextstep, force_LJ)

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



    #return simLog
    return poslog
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


    
#end # module