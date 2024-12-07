########## THIS FILE WILL NOT USE NAIVEDYNAMICS DEV ENVIRONMENT. IT WILL USE THE GITHUB VERSION
using Revise
using Test

#using BenchmarkTools
#using CSV
#using StaticArrays
using NaiveDynamics
#import GLMakie
#Base.get_extension(NaiveDynamics, NaiveMakie)
#using GLMakie
#using oneAPI
#using oneAPI.oneL0


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



myCollector1 = GenericRandomCollector(; floattype=Float32,
                                    objectnumber=8,
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
myCollection1 = collect_objects(myCollector1)

bvhspec = SpheresBVHSpecs(; floattype=Float32, 
                            interaction_distance=0.1, 
                            atoms_count=length(myCollection1.position), 
                            bins_count=length(myCollection1.position) 
)
build_bvh(myCollection1.position, bvhspec, myCollector1 )


#myCollector = GenericRandomCollector(Float32, 40, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -0.02, 0.02, 0.001)
#myCollection = collect_objects(myCollector)
#mySpec = GenericSpec{Int64}(4000, 1, 1, 10)
#logpos2 = simulate_dumloop!(myCollection, mySpec, myCollector)
#@btime logpos2 = simulate!($myCollection, $mySpec, $myCollector)
