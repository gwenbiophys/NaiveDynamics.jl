# NaiveDynamics.jl
[![Build status (Github Actions)](https://github.com/gwenbiophys/NaiveDynamics.jl/workflows/CI/badge.svg)](https://github.com/gwenbiophys/NaiveDynamics.jl/actions)
[![codecov](https://codecov.io/github/gwenbiophys/NaiveDynamics.jl/graph/badge.svg?token=MMODZ51EE5)](https://codecov.io/github/gwenbiophys/NaiveDynamics.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://gwenbiophys.github.io/NaiveDynamics.jl/dev)


A very Naive Project to follow and study the design of molecular dynamics simulations and high octane algorithms.

As of 13 March, 2025:
* unitless velocity verlet simulation 
* no-identity point-particles
* partially complete coulomb and Lennard-Jones force models
* an eerily unhelpful velocity rescaler
* ability to video record simulation
* bvh traversal for neighbor search (CPU parallel) -- it actually works now !
