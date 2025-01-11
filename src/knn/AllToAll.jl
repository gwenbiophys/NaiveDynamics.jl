export
    update_pairslist!,
    unique_pairs,
    threshold_pairs,
    new_threshold_pairs

@inline function pairslist_interior(each, a::Vec3D{T}, list) where T

    i = list[each][1]
    j = list[each][2]
    dxyz = SVector(a[i] .- a[j])
    d2 = sqrt( sum(dxyz .^ 2))



    list[each] = tuple(i, j, dxyz, d2)



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
    # result = tuple(i, j, dx, dy, dz, d2)


    # list[each] = result


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
    list = [tuple(i, j, SVector{3, T}(0.0, 0.0, 0.0), a[1][1]) for i in 1:length(a)-1 for j in i+1:length(a)]
    #list =  [tuple(i, j, a[1][1]) for i in 1:length(a)-1 for j in i+1:length(a)]
    #fillcombinations


    update_pairslist!(a, list)


    return list
end

#so. how can we improve this, as it takes 84% of execution time
# I don't think Julia likes it lol
@inline function threshold_pairs(list, threshold::T) where T
    
    return [list[i] for i in eachindex(list) if list[i][4] ≤ threshold] #TODO this is here hrmm

end

#TODO this is MUCH more mrmory efficient, but is also slightly universally slower from 10 to 5000
# with @time, usually this loses, but almost always loses to @btime
@inline function new_threshold_pairs(list, threshold::T) where T
    counter = 0
    for each in eachindex(list)
        if list[each][3] ≤ threshold
            counter += 1
        end
    end
    thresh = [tuple(0, 0, 0.0f0) for i in 1:counter]#Vector{Tuple{Int64, Int64, T}}(undef, counter)
    #resize!(thresh, counter)
     #hrm
    for each in eachindex(thresh)
        if list[each][3] ≤ threshold
            thresh[each] = list[each]
        end
    end
    return thresh#[list[i] for i in eachindex(list) if list[i][3] ≤ threshold] 

end
