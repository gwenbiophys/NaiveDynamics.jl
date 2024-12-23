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
#using CellListMap
#using NearestNeighbors
using SortingAlgorithms

using Aqua #is this correct? should probably go in test files!

#using KernelAbstractions
#using AMDGPU
#using CUDA
#using GLMakie
#using InteractiveUtils


include("MDInput.jl")

include("BaseExtensions.jl")
#include("knn/NeighborSearch.jl")
include("knn/ProkoLebrun.jl")
include("PkgExtensions.jl")
include("knn/AllToAll.jl")
include("Forces.jl")
include("Simulator.jl")


end # module NaiveDynamics
