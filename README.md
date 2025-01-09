# NaiveDynamics.jl
[![Build status (Github Actions)](https://github.com/gwenbiophys/NaiveDynamics.jl/workflows/CI/badge.svg)](https://github.com/gwenbiophys/NaiveDynamics.jl/actions)
[![codecov](https://codecov.io/github/gwenbiophys/NaiveDynamics.jl/graph/badge.svg?token=MMODZ51EE5)](https://codecov.io/github/gwenbiophys/NaiveDynamics.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://gwenbiophys.github.io/NaiveDynamics.jl/dev)


A very Naive Project to follow and study the design of molecular dynamics simulation and high octane algorithms.

As of 9 January, 2025:
* unitless velocity verlet simulation 
* no-identity point-particles
* partially complete coulomb and Lennard-Jones force models
* an eerily unhelpful velocity rescaler
* ability video record simulation
* bvh traversal for neighbor search (CPU parallel)


<video controls src="https://github.com/gwenbiophys/NaiveDynamics.jl/blob/main/data/newhope.mp4" title="SampleSim"></video>