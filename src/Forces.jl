export
    force_lennardjones!,
    force_coulomb!,
    sum_forces!

Base.@propagate_inbounds function lennardjones_interior(each, eps, σ, force::Vec3D{T}, pairslist) where T
    #d2 = pairslist[each][3]
    # TODO
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
    



    # incorrect, overall force is being misapplied to each component
    #force[i] .+= (24*eps ./ d ) .* ((2*σ ./ d).^12 .- (σ ./ d).^6)
    #force[j] .-= force[i]
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

    Threads.@threads for each in eachindex(pairslist)
        @inbounds lennardjones_interior(each, eps, σ, force, pairslist)
        
    end


    #return force
end
Base.@propagate_inbounds function coulomb_interior(each, k, force::Vec3D{T}, pairslist, charge) where T
    i = pairslist[each][1]
    j = pairslist[each][2]
    dx = pairslist[each][3]
    dy = pairslist[each][4]
    dz = pairslist[each][5]

    #d = pairslist[each][6]

    force[i][1] += k * charge[i] * charge[j] / (dx^2)
    force[i][2] += k * charge[i] * charge[j] / (dy^2)
    force[i][3] += k * charge[i] * charge[j] / (dz^2)
    force[j][1] -= force[i][1]
    force[j][2] -= force[i][2]
    force[j][3] -= force[i][3]
end


function force_coulomb!(force::Vec3D{T}, pairslist, charge) where T
    k = 1 # coulomb constant
    for each in eachindex(force)
        map!(x->x, force[each], MVector{3, T}(0.0, 0.0, 0.0))
    end

    Threads.@threads for each in eachindex(pairslist)
        @inbounds coulomb_interior(each, k, force, pairslist, charge)

    end
end

function sum_forces!(force, force1, force2)

    for each in eachindex(force)
        force[each] .= force1[each] .+ force2[each] 
    end

    return force
end