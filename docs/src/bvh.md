# Implementing Bounding Volume Hierarchies in NaiveDynamics.jl
I am working towards a stackless BVH traversal algorithm for accelerating neighbor search in molecular dynamics simulations. This effort is attempting to implement A. Prokopenko and Lebrun-Grandie's method [^1].



## Current Roadmap
- [x] Reimplement in Julia
- [] Optimize performance and improve code readability/structure
- [] Benchmark and report my findings here [^3]
- [] Build a comparable testing situation between ArborX and NaiveDynamics 

## Motive
1. Because performance
    With modern graphics processors adding ray tracing hardware for realistic rendering, we have a *small* opportunity to plug these hardware blocks into molecular dynamics simulations. For instance, Intel's 'Xe2' graphics architecture can run 18 box tests per clock per ray-tracing unit [^2], and do so indepentently of the vector processors. That could be a wonderful tool.
    


2. Because traditional Verlet lists are not appropriate for all simulation types

3. because it is cool! What better way could there be to be introduced to GPU programming than trying on something that might not work?



## Journey towards testing









## Footnotes

[^1]: [Revising Apetrei’s bounding volume hierarchy construction algorithm to allow stackless traversal] (https://info.ornl.gov/sites/publications/Files/Pub208673.pdf)
[^2]: [Chips and Cheese: Lunar Lake iGPU] (https://chipsandcheese.com/p/lunar-lakes-igpu-debut-of-intels)
[^3]: Testing was conducted on my machine( i5-9300H, stock CachyOS installed in November, 4 Threads enabled) using Base.@showtime and BenchmarkTools.@btime