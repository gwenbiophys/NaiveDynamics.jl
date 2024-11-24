using Revise
using Test

#using BenchmarkTools
#using CSV
#using StaticArrays
using NaiveDynamics
import GLMakie

@testset begin
    myCollector = GenericRandomCollector(; floattype=Float32,
                                        objectnumber=50,
                                        minDim=tuple(-1.0, -1.0, -1.0),
                                        maxDim=tuple(1.0, 1.0, 1.0),
                                        temperature=0.01,
                                        randomvelocity=false,
                                        minmass=1.0,
                                        maxmass=5.0,
                                        minimumdistance=0.001,
                                        mincharge=-1f-9,
                                        maxcharge=1f-9
                                        )

    myCollection = collect_objects(myCollector)
    #mySpec = GenericSpec{Int64, Float32}(50, 1, 1, 10, 1)
    mySpec = GenericSpec(; inttype=Int64,
                        floattype=Float32,
                        duration=10,
                        stepwidth=1,
                        currentstep=1,
                        logLength=10,
                        vDamp=1)
    logpos = simulate!(myCollection, mySpec, myCollector)
    #@profview simulate!(myCollection, mySpec, myCollector)
    #@btime logpos2 = simulate!($myCollection, $mySpec, $myCollector)
    direc = "/home/gwenk/Coding/Julia/NaiveDynamics.jl/data/iWant.mp4"
    record_video(direc, logpos, myCollector; frameinterval = 1)
end