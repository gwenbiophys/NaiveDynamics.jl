export
    update_pairslist!,
    unique_pairs,
    threshold_pairs

@inline function pairslist_interior(each, a::Vec3D{T}, list) where T

    i = list[each][1]
    j = list[each][2]
    d2 = sqrt( sum((a[i] .- a[j]) .^ 2))
    # xi = a[i][1]
    # xj = a[j][1]

    # yi = a[i][2]
    # yj = a[j][2]

    # zi = a[i][3]
    # zj = a[j][3]

    # dx = xi - xj
    # dy = yi - yj
    # dz = zi - zj
    # d2 = sqrt(dx^2 + dy^2 + dz^2)  
    #result = tuple(i, j, dx, dy, dz, d2)
    list[each] = tuple(i, j, d2)

    #list[each] = result


end
function update_pairslist!(a::Vec3D{T}, list) where T

    #TODO this if should not be evaluated here, should be evaluate4d ahead of time
    if length(a) == 2

        dx = a[1][1] - a[2][1]
        dy = a[1][2] - a[2][2] 
        dz = a[1][3] - a[2][3]
        d2 = sqrt(dx^2 + dy^2 + dz^2) 

        list[1] = tuple(1, 2, d2)

    else

        Threads.@threads for each in eachindex(list)
            pairslist_interior(each, a, list)
        end
    end

    return list
end

function unique_pairs(a::Vec3D{T}) where T
    #list_length = convert(Int64, (length(a)-1) * length(a) / 2)

#TODO could we do another thread chunking scheme where we initialize and combine them at the end?
#or is current scheme ideal, where we initialize dummy values and then parallelwise fill them in later?
    #list = [tuple(i, j, a[1][1], a[1][1], a[1][1], a[1][1]) for i in 1:length(a)-1 for j in i+1:length(a)]
    list =  [tuple(i, j, a[1][1]) for i in 1:length(a)-1 for j in i+1:length(a)]
    #fillcombinations


    update_pairslist!(a, list)


    return list
end

#i don't think this is worth parallelizing
function threshold_pairs(list, threshold::T) where T
    
    return [list[i] for i in eachindex(list) if list[i][3] ≤ threshold] #TODO this is here hrmm

end
