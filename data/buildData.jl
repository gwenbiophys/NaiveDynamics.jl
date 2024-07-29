#using Revise
#using BenchmarkTools
#using CSV
#using StaticArrays
#using GLMakie


#using ProfileView #doesnt work in VSCode




#myCollector = GenericRandomCollector{Float32}(40, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -0.000002, 0.000002, 0.001)
#myCollection = collect_objects(myCollector)
#mySpec = GenericSpec{Int64}(800000, 1, 1, 10)
# memory limit starts getting feisty at 4 million steps
#@btime simulate!($largeSimulation, $largeCollector)
#@btime logpos = simulate_bravado!($myCollection, $mySpec, $myCollector)
#@code_warntype simulate_bravado!(myCollection, mySpec, myCollector)
#@profview logpos = simulate!(myCollection, mySpec, myCollector)
#myCollector = GenericRandomCollector{Float32}(50, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 0.0001, false, 1.0, 5.0, 0.001, -5.0, 5.0)
myCollector = GenericRandomCollector(; floattype=Float32,
                                    objectnumber=50,
                                    min_xDim=-1.0,
                                    min_yDim=-1.0,
                                    min_zDim=-1.0,
                                    max_xDim=1.0,
                                    max_yDim=1.0,
                                    max_zDim=1.0,
                                    temperature=0.0001,
                                    randomvelocity=false,
                                    minmass=1.0,
                                    maxmass=5.0,
                                    minimumdistance=0.001,
                                    mincharge=-5.0,
                                    maxcharge=5.0
                                    )
myCollection = collect_objects(myCollector)
#mySpec = GenericSpec{Int64, Float32}(50, 1, 1, 10, 1)
mySpec = GenericSpec(; inttype=Int64,
                    floattype=Float32,
                    duration=50,
                    stepwidth=1,
                    currentstep=1,
                    logLength=10,
                    vDamp=1)
@time logpos = simulate!(myCollection, mySpec, myCollector)
#@btime logpos2 = simulate!($myCollection, $mySpec, $myCollector)

@time record_video("C:/Users/kucer/Desktop/julia/NaiveDynamics.jl/data/newhope.mp4", logpos, myCollector; frameinterval = 1 )





#myCollector = GenericRandomCollector(Float32, 40, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -0.02, 0.02, 0.001)
#myCollection = collect_objects(myCollector)
#mySpec = GenericSpec{Int64}(4000, 1, 1, 10)
#logpos2 = simulate_dumloop!(myCollection, mySpec, myCollector)
#@btime logpos2 = simulate!($myCollection, $mySpec, $myCollector)
