export
    force_lennardjones!,
    force_coulomb!,
    sum_forces!

Base.@propagate_inbounds function lennardjones_interior(each, eps, σ, force::Vec3D{T}, pairslist) where T
    #d2 = pairslist[each][3]
    # TODO
    i = pairslist[each][1]
    j = pairslist[each][2]
    dxyz = pairslist[each][3]
    force[i] .+= (24*eps ./ dxyz ) .* ((2*σ ./ dxyz) .^ T(12.0) .- (σ ./ dxyz) .^ T(6.0))

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
        @inbounds lennardjones_interior(each, eps, σ, force, pairslist)
        
    end


    #return force
end
Base.@propagate_inbounds function coulomb_interior(each, k, force::Vec3D{T}, pairslist, charge) where T
    i = pairslist[each][1]
    j = pairslist[each][2]
    dxyz = pairslist[each][3]
    force[i] .+= k * charge[i] * charge[j] ./ (dxyz .^ 2)
    force[j] .-= force[i]

end


function force_coulomb!(force::Vec3D{T}, pairslist, charge) where T
    k = 1 # coulomb constant
    for each in eachindex(force)
        map!(x->x, force[each], MVector{3, T}(0.0, 0.0, 0.0))
    end

    for each in eachindex(pairslist)
        @inbounds coulomb_interior(each, k, force, pairslist, charge)

    end
end

function sum_forces!(force, force1, force2)

    for each in eachindex(force)
        force[each] .= force1[each] .+ force2[each] 
    end

    return force
end