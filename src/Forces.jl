export
    update_pairlist!,
    unique_pairs,
    threshold_pairs,
    force_lennardjones!,
    force_coulomb!,
    sum_forces!


Base.@propagate_inbounds function pairslist_interior(each, a::Vec3D{T}, list) where T

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
function update_pairslist!(a::Vec3D{T}, list) where T

    
    if length(a) == 2

        dx = a[1][1] - a[2][1]
        dy = a[1][2] - a[2][2] 
        dz = a[1][3] - a[2][3]
        d2 = sqrt(dx^2 + dy^2 + dz^2) 

        list[1] = tuple(1, 2, dx, dy, dz, d2)

    else


        Threads.@threads for each in eachindex(list)
            @inbounds pairslist_interior(each, a, list)
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

Base.@propagate_inbounds function lennardjones_interior(each, eps, σ, force::Vec3D{T}, pairslist) where T
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