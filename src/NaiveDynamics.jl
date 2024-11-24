module NaiveDynamics

using BenchmarkTools
using UUIDs
using CSV
using DataFrames
using NamedArrays
using StaticArrays
using Distributions
using Revise
using StructArrays
using CellListMap
using NearestNeighbors

using Aqua #is this correct? should probably go in test files!

#using KernelAbstractions
#using AMDGPU
#using CUDA
#using GLMakie


include("MDInput.jl")
include("Simulator.jl")
include("BaseExtensions.jl")
include("NeighborSearch.jl")
include("PkgExtensions.jl")
include("Forces.jl")


end # module NaiveDynamics
