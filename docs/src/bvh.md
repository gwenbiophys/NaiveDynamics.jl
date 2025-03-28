# Implementing Bounding Volume Hierarchies in NaiveDynamics.jl
I am working towards a stackless BVH traversal algorithm for accelerating neighbor search in molecular dynamics simulations. This effort is attempting to implement A. Prokopenko and Lebrun-Grandie's method [^1].



## Current Roadmap
* [x] Reimplement in Julia
* [] Optimize performance and improve code readability/structure
* [] Benchmark and report my findings here [^2]
* [] Build a comparable testing situation between ArborX and NaiveDynamics 

## Motive
1. Because performance
With modern graphics processors adding ray tracing hardware for realistic rendering, we have a *small and theoretical* opportunity to plug these hardware blocks into molecular dynamics simulations. For instance, Intel's 'Xe2' graphics architecture can run 18 box tests per clock per ray-tracing unit [^3], and do so indepentently of the vector processors. That could be a wonderful tool in a world where my own CPU execution of bvh building, traversal, and a dummy simulation sees box testing taking 60% of total execution time.
    


2. Because traditional Verlet lists are not appropriate for all simulation types.

3. Because it is cool! What better way could there be to be introduced to programming than trying something I have no idea how to do!



## Journey towards testing









## Footnotes

[^1]: [Revising Apetrei’s bounding volume hierarchy construction algorithm to allow stackless traversal] (https://info.ornl.gov/sites/publications/Files/Pub208673.pdf)
[^2]: Testing was conducted on my machine( i5-9300H, stock CachyOS installed in November, 8 Threads enabled)
[^3]: [Chips and Cheese: Lunar Lake iGPU] (https://chipsandcheese.com/p/lunar-lakes-igpu-debut-of-intels)
