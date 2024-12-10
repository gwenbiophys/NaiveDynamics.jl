##### data and tester funcs

function self_traverse(keys, neighbors, position, spec::SpheresBVHSpecs{T, K}) where {T, K}
    currentKey = 1 + spec.leaves_count
    selfVisitedLeaves = [false for each in eachindex(position)]

    while true



        if currentKey == 0 # is sentinel node
            break
        end
    end

    if sum(selfVisitedLeaves) == spec.leaves_count
        return true
    else
        return false
    end


end
myCollector8 = GenericRandomCollector(; floattype=Float32,
                                    objectnumber=8,
                                    minDim=tuple(0.0, 0.0, 0.0),
                                    maxDim=tuple(1.0, 1.0, 1.0),
                                    temperature=0.01,
                                    randomvelocity=false,
                                    minmass=1.0,
                                    maxmass=5.0,
                                    minimumdistance=0.001,
                                    mincharge=-1f-9,
                                    maxcharge=1f-9
)
position8 =[MVector{3, Float32}(0.1, 0.1, 0.1), MVector{3, Float32}(0.2, 0.2, 0.2), 
            MVector{3, Float32}(0.346, 0.98, 0.12), MVector{3, Float32}(0.01, 0.76, 0.99), 
            MVector{3, Float32}(0.1111, 0.4, 0.31), MVector{3, Float32}(0.234, 0.29, 0.2), 
            MVector{3, Float32}(0.11346, 0.918, 0.1276), MVector{3, Float32}(0.061, 0.76, 0.989)
]
bvhspec8 = SpheresBVHSpecs(; floattype=Float32, 
                            interaction_distance=0.1, 
                            leaves_count=length(position8) 
)
keys = build_bvh(position8, bvhspec8, myCollector8)


##### functionality tests
@testset "traversability" begin

    @test keys[8].skip == 0
    for each in 1:7
        @test keys[each].skip != 0
    end
end

@testset "root boundaries include all atoms" begin
    for each in 1:7
        for i in eachindex(position8[1])
            @test keys[8].min[i] < position8[each][i] < keys[8].max[i]
        end
    end
end

# @testset "morton sorting" begin

# end


##### performance tests(?)
