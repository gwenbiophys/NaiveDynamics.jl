# Developer Diary

## 1. Let's document - 27 April
Today has been an interesting day. I have been working on my baseline code but really hardening my documentation so that I have a well collected space to start hucking my thoughts. 

For instance, I added CSV and NamedArrays .jl to prepare for testing on whether the 
```julia
mutable struct GenericObjectCollection <: ObjectCollection 
    name::AbstractArray{String, 1}

    position::AbstractArray{AbstractFloat, 3}
    velocity::AbstractArray{AbstractFloat, 3}

    uniqueID::AbstractArray{UUID,1}

end
```
makes any sense, or if I should place all this information into a single array. And then I could test how this version of GenericObjectCollection scales, comparing it against a NamedArray convention. The primary point is to minimize processing time on putting these arrays together, while keeping my function accesses to data meaningful. I do not want the following:

```julia
Collection = [AbstractArray{String, 1}, AbstractArray{AbstractFloat, 3}, AbstractArray{AbstractFloat, 3}]
function do_something(Collection)
    return Collection[1] - Collection[2]
end
```
in which function writing depends on making sure I have the right index of my Collection, so I don't do something stupid, like add a force term to a velocity term. I admit, something like ``` position = Collection[3]``` is reallly easy, but we are here for overengineered solutions. ;)


But then my julia package would not precompile! Because I added new packages without updating the project.toml, so I had to relearn how to add new dependencies to a package.



It has been an hour and I still haven't figured out how to automatically get docstrings for functions passed into an index/api page, like they have it over at Molly.jl. Oh well, it is likely far better to get the gh-pages version of the documentation working.

## 2. 3 May - packaging packages locally

in order to point a test file to my development package and skipall the githubbing and comppiling, i must
```julia
pgk() dev ./NaiveDynamics.jl
using Revise
using NaiveDynamics
```
in that specific order. That way a change can be made in the development-source, and immediately referenced in the test/user file. Thank you very much, Revise.jl/stable. To hell with you, Copilot.



### Let's talk about DataFrame-ification
In the commits last weekend I was very focused on modifiying this package's structure so that a user, regardless of usecase, could design functions around their data in the ```ObjectCollection``` type WITHOUT relying on the order of set items in the type. In the current ```collect_objects()``` function, I have this nonsense:
```julia 
    simCollection = GenericObjectCollection(
        fill(step_n, objectcount),
        fill("duck", objectcount),
        [1:objectcount;],
        [SizedVector{3, Float64}(rand(positionRange, 3)) for each in 1:objectcount],
        [SizedVector{3, Float64}(rand(velocityRange, 3)) for each in 1:objectcount],
        [SizedVector{3, Float64}(zeros(Float64, 3)) for each in 1:objectcount],
        )
```
when I would prefer a semantics based method that does not depend on the order of inputs. After a lot of testing I just tossed in a big dataframe and a little dataframe. I ran a 1 step simulation of 100 000 atoms, and execution time took at least 3.7 seconds and processed through a humonculous 2GiB of data. I tweaked it a little, and then it wouldnt run anymore at all for OutOfMemory errors. So it became my new task to performance test again and again on all the different methods and packages I could find for multiplying two arrays by broadcasting. We got to a very good place with generating position vectors:
```julia
position = [MVector{3, Float64}(rand(posRange, 3)) for each in 1:5]
```
and with performing the multiplication, as base.broadcast has no idea how to solve a vector of mutable vectors:
```julia
function posVel_multiply!(position, velocity)
    for i in eachindex(position)
        position[i] .*= velocity[i]
    end
    return position
end
```
We had a path to de-dataframe-ification and a massive performance improvement. But we still had not solved the convenience problem, but ya know what? I'm good with inconvenient. I have spent a long time trying to solve a problem that isn't a problem. I can winback my semantics with statements at the beginning of a function, like in simulate where we unwind the datatype hierarchy to pull the position array of arrays out of the "system" type. It's a less beautiful solution maybe, but no less effective.

And how much better is life? Sans parallelism and multithreading--just from algorithmic and datatype-selection improvements-- we need to run the 100k simulation out to ~60 steps for it to take as long as the DataFrame's version. And we need to run at least 130 steps to generate as much information as the DataFrame's 1 step. Not that I blame the package's authors, it was a very silly idea to turn a highly versatile data processing tool into the DNA of the organism that forms this Naive package. 

There are still algorithmic improvements to be made. I wanted to use a vector of mutable vectors, like this:
```julia
using StaticArrays
using Distributions

minposition = -5.0
maxposition = 5.0
posRange = Uniform(minposition, maxposition)
velRange = Uniform(minposition, maxposition)
position = [MVector{3, Float64}(rand(posRange, 3)) for each in 1:5]
velocity = [MVector{3, Float64}(rand(velRange, 3)) for each in 1:5]


function posVel_multiply!(position, velocity)
    for i in eachindex(position)
        position[i] .*= velocity[i]
    end
    return position
end

posVel_multiply!(position, velocity)
```
but the posVel_multiply! does not work in the simulate!() function, and I have absolutely no clue why. It looks to be an active issue on github, and substituting the mutablevector with a sized vector solves the issue. But at low atoms counts, performance falls back by about 20%, and it balloons as count increases. At 100k atoms, it appears about 10x slower than a mutable vector, at about 3ms for the broadcasted multiplication alone. In context, the sized vector can only multiply the position and velocity vectors in about 50 ms per step, and this is an engine with zero interactions. 

However the logging of the simulation actually takes longer than everything else under ``` simulate!() ```, so there are several angles of optimization before we even consider threading or SIMD operations that the compiler doesn't already perform.

But for now we have an angle of approach a velocity verlet stepper that will help bring the time per step back up to >1 second

### Towards a VerletVelocity Lennard-Jones fluid
The present implementation of ```unique_pairlist``` scales at least with BigO(n^2), drawing out execution time for 1 step at 100k particles to over 20 minutes.
So the CellListMap method will be used to calculate neighbors instead until I either learn about how it works or come up with some Naive mapping routine. I wonder if there is some shenanigans to pull around by drawing a weirdly shaped matrix, and in so instantiating it, we will determine all of the unique pairs in the our dataset, then we just have to redraw the matrix into a very long vector. However, maybe the real magic of the CellListMap method is the ability to avoid first defining all unique pairs and then asking which ones are within a given distance. And instead it somehow directly creats a list of unique pairs using the distance as the 'first' metric, rather than the last metric. Whatever it is I do want to continue investigating this and create a real Naive method, just to explore the algorithms and the sense of going from
```julia
    for i in eachindex(position)
        for j in eachindex(position) 
            if j > i
                push!(list, SVector(position[i],position[j]))
            end
        end
    end
```
to something not only more useful, but requiring different coding techniques.

I haven't an idea how these would even work, so for now let's just make CellListMap operable.So we are not quite making CellListMaps, but neighbor lists instead. Not as fancy but still excellent.


If the ```force_currentstep``` is the sum of each force acting on/within the system, must I necessarily allocate and deallocate an array for each force and then broadcast-add them to ```force_currentstep```, or could I instead perform each force-calculation using ```force_currentstep``` as an argument and returning it?  I believe the second method is more memory efficient but less parallelizable, when the first method is parallelizable until my device runs out of memory.
Naively, doing work sequentially would look either like a series of equal signs
```julia
force_currentstep = LennJones_force!(force_currentstep)
force_currentstep = coulomb_force!(force_currentstep)
etc
```
or as a nested set of functions
```julia
force_currentstep = coulomb_force!(LennJones_force!(force_currentstep))
```
For the sake of difference, let's use the nested set one :) . Here's an ambition: a simulation-time optimization-suite that detects pending out of memory error/is system aware and decides if we have to focus on memory management or we can go wild in parallel

## 3. 12 May - curious learning moment
So today I have overhauled how positions are generated to allow the user to define the size of the box with different sizes in each dimension. This is in preparation for defining a reflective box boundary later on. And we are in a strange situation, when calling the code 
```julia
pairslist = InPlaceNeighborList(x=position, cutoff=0.1, parallel=false)
```
in my simulate! function, CellListMap.jl checks the position values' maximum and minimum according to their datatype. This check fails because for whatever reason, the number type in the position vector of vectors is an AbstractFloat, despite verification demonstrating that there are real values here. What is even stranger is that I have asked the datatype of my position numbers immediately after creation, and they are Float64. I even tried changing the definition in the constructor GenericRandomCollector from AbstractFloat's to Float64's. 

Here is the list of things that separate the position vector from it's point of creation to the point of being bassed on to INPL: 
1. The vecOfMVecs is created under ```generate_positions(Collector)``` and is correctly typed as Float64.
2. The vecOfMVecs is retured from ```generate_positions(Collector)``` to the positions spot where simCollection is initialized in ```collect_objects(Collector)``` and is incorrectly typed as AbstractFloat. It appears that the initialization of a GenericObjectCollection datatype is corrupting the Float64 to AbstractFloat. This seems as unexpected behavior to me, why would it do this? Redefining the GenericObjectCollection fixes the issue. But why would a concrete type ever be reverted back to an abstract type?

Because Julia is a bit sneaky with type promotion. I misremembered that Julia never promoted types from concrete to abstract. A way of solving this issue, I have learned, is to introduce parametric typing into my code. And this is very exciting! I didn't quite understand any applications of parametric typing (mostly for lack of effort) and here is one just for me! Parametric typing to prevent promotion, moreover, to perform type inferencing based only on the user input.  

### Reworking broadcasting and operators in velocity verlet method
Because we have to work with neighbor listing and mapping to calculate our forces, it's a bit harder to ```simulate!()``` by having a nested loop going over every atom. So instead we have to rewrite the calculation of each line somehow. I tried various methods
1. Broadcast-broadcast by using a unique syntax
``` ⊕(a, b) = a .+ b``` and then calling ```.⊕``` to reach the data points within our VecOfVecs
2. Mapping a broadcasted-operator to our VecOfVecs
```map!(.*, $a,$b)``` and ```map!(./, $a,$b, $a)```
3. Dumloop over each item of the outer-vector, and then broadcast an operation to each item
```julia
function dumloop_multiply!(d, e)
    @inbounds for i in eachindex(d)
        d[i] .*= e[i]
    end
end
```
The I tried to make a version with an *operator* argument, but to no success.
```
function dumloop_multiply!(operator d, e)           function dumloop_multiply!(operator d, e)
    @inbounds for i in eachindex(d)                     @inbounds for i in eachindex(d) 
        d[i] .operator= e[i]                                d[i] = operator(d[i], e[i])
    end                                                 end
end                                                  end
```
The goal is to extend the behaviors of arithmetic operators and broadcasting to work directly on Vectors of Vectors. However, maybe I can avoid this nonsense by working on how the force is calculated. And I have done that. 

### JL force
And I quickly ran into a new problem: The neighborlist method in CellListMap.jl returns a vector of tuples ```(indexA, indexB, distance_as_the_fish_swims)```, when my algorithms demand vectorized data. I need a neighborlist algorithm to return me a tuple of ```indexA, indexB, MVector(distanceX, distanceY, distanceZ), distanceLine``` in order to correctly calculate the component forces and the overall energy of these interactions. Or maybe I can code that myself hahah!

This would involve using the CellListMap interface, but I have a streaming method in which each atom is given a block of calculations to run through, so parallelism would primarily come from the top forloop by calculating multiple atoms at once in more or less the same procedure. So I **really** need to make a naive neighbor list with all of the data I need, or otherwise figure out how to pull that information out of CellListMap, but that place is designed to scare people.

A further quandary arises from the highly differing behaviors of neighborlist() and the ```InPlaceNeighborList()```, ```neighborlist!()```, and ```update!()``` which is the latter does not return the neighbor list vector of tuples, but a complicated data structure series. It is possible that this series contains my data of xyz distances, but sheesh. It may be happier if I provided box-boundaries, as the compute engine takes a boundary that is either user determined or limited by the values expressable in the data-type of the coordinate (Float64 makes the InPlaceNeighborList object fill up the entire REPL with almost only zeros, we learned from our type promotion earlier). Let's see how slow the sim engine is now, and maybe we might consider making a switch to NearestNeighbors.jl


### the simulation BREAKS
Testing the package as it presently stands, it works fine for 100 atoms, but the upperbound number is falling as I examine it by hitting the mid point between the working number and the failing number. Currently, the simulation crashes out at or below 600 atoms. This is a method error, ```MethodError: objects of type StaticArraysCore.MVector{3, Float64} are not callable``` that points to these two lines in Simulator.jl:
```julia
d = position[i] .- position[j]
force[i] .+= (24*eps ./ d )((2*σ ./ d).^12 - (σ ./ d).^6)
```
WAIT, the fact that some runs show problems and others do not indicate that their atoms are not within the cutoff. It's a gentle bug coming from my work! Very gentle, here's the fixed function: 
```julia 
force[i] .+= (24*eps ./ d )((2*σ ./ d).^12 - (σ ./ d).^6)
```


Onto the next breakage ```InexactError: trunc(Int64, 3.884345474457886e38)```. 
A more severe one, I fear. It points to the ```neighborlist()``` function call. We are going to update CellListMap and pray. A new release came about 3 days ago sooo we should be good, right? *Oh dear.* Going from 1000 to 10 000 particles seems to drastically increase the size of the error here, from about 10^20 to 10^30 through and beyond 10^50. So we have to figure out what value is being created that is so extremely stratospheric, and it is created seemingly each time enough particles are created to exist in less than the cut off. Let's increase cut off and reduce atoms. Now we are getting different errors regardless of using a function or sitting in global scope, outofmemory, invalid array size, invalid array dimensions, and sometimes the truncation error. Each of these errors occur in the neighborlist run. For the truncation error, the problematic calculation is Line 213 of Box.jl. 
```_nc = floor.(Int, (xmax .- xmin) / (cutoff / lcell)) 
``` 
NearestNeighbors.jl hass it's own problem that the most recent issue shows a method error that prevents the use of Vector(MVector()) so that may not be a solution.

Ah but wait, I did not correctly construct my force updates. Now it is working! I am uncertain how I was obliterating the neighborlist by not calculating forces correctly, but it is a strange world. Now it is not working anymore, after I increased the stepcount. Failing in the neighbor list for multiple different reasons. There are 3 erros, out of memory, invalid arraysize, invalid array dimenions, and they all are genereated in line 482 of CellLists.jl
```julia
    cl = CellList{N,T}(
        n_real_particles=length(x),
        number_of_cells=prod(box.nc),
    )
```
It is almost as though the REPL gets 'gunked' up from multiple rounds of ```simulate!()``` Or that perhaps at times two particles land on top of each other or there is some other issue. The gunking theory fell out because no change resulted in a successful run.

We have 2 issue classes, truncation from Box.jl 213 and CellLists.jl 482.

Wait a minute, perhaps the particles are falling in too close and the forces are spiraling out of control. Changing epsilon by 1000x either way seemed to have little effect. While reducing sigma by 1000x seemed to improve success and increasing seemed to mkae failure fairly common. However, success is not guaranteed within 100 or 1000 runs. Hopefully with a pruning routine, things will turn out just fine. An additional helper is increasing the cutoff. Perhaps a class of errors, or maybe jsut a specifc error or a few can be eliminated with a high enough cut-off or dense enough conditions to guarantee that the neighbor list is not empty.

19 May
Even with the changes, the naive pruning and the like, we are still getting NotANumber answers, which may indicated an enduring problem with parameterization in the neighbor cutoff and in the time width and in the particle density and etc. But at least the simualtion runs without dying now, if it runs inconsistently due to the use of in-Verlet try catch blocking. 

A fairly curious discussion comes around strategies for conserving energy. At present, our work is purely additive. Seemingly after the first step, all particles will be slammed into the walls of the simulation box.


## 4. 24 May - Death to allocations

So my issue with CellListMap is a limitation of the method as a whole, whose memory scales with the volume of space being considered, or something like that, and my highly disperse simulations seem to force the algorithm to consider the whole numeric width of Float64 for calculating the neighbor list. Thanks to the university of michigan for their informative article.


Before I worry about pulling in NearestNeighbors.jl or seeing how long it would take to make a naive BVH neighbor finder, I am allocation hunting in my simulation to increase speed. I was primarily worried that ```force_currentstep = force_nextstep``` would de/allocate, as I often run into the problem that ```a = 5, b = a, b += 1, but a !== 6```, but the compiler knows what it's doing. We instead run into a much more interesting problem:

To perform Velocity-Verlet integration, we use broadcasting to calculate the xyz-components of position for each particle:
```julia
position[i] = position[i] .+ (velocity[i] .* stepwidth) .+ (force_currentstep[i] ./ mass[i] .* stepwidth^2/2)
```
Benchmarking showed an allocation was made, so I set to improve it:
```julia
positionJ[i] .+= ((velocity[i] .* stepwidth) .+ (force_currentstep[i] ./ mass[i] .* stepwidth^2/2))
```

The two methods should obviously (to me) produce the same output, and are seemingly just format differences. However when calculating positionJ, we increase speed, 16 -> 10 ns, and eliminate allocation, 32 bytes -> 0. That alone, positionJ seems to be the idea choice. But why the difference?

More importantly, what about correctness? positionJ and position produce different results about 5% of the time. It turns out positionJ's method is internally consistent, while position's is not internally consistent. And that is wild.

We find internal consistency by either removing the velocity-stepwidth term, or by using positionJ's method. Removing the addition of ```(force_currentstep[i] ./ mass[i] .* stepwidth^2/2)``` to position does not affect internal consistency. Furthermore, switching the order of operations
```julia
position[i] = position[i]  .+ (force_currentstep[i] ./ mass[i] .* stepwidth^2/2) .+ (velocity[i] .* stepwidth)
positionK[i] = positionK[i] .+ (velocity[i] .* stepwidth) .+ (force_currentstep[i] ./ mass[i] .* stepwidth^2/2)
```
also has no affect, position and positionK's methods are still inconsistent with each other.

So, I just made my code faster* and more precise by changing how it was spelled. What a weird world. *A savings of 6 ns per calculation doesn't register at the ms level in my short-small test, but memory usage and allocations are reduced by 6.7 and 7.7 percent, respectively. 

The next line of question, why is simulate!() (with neighborlisting and force-calculating deactivated) so slow and allocating so often? Internally, the function should allocate only a few times, with each action taking fewer than 50 ns, and most actions taking fewer than 20 ns. Changing the index scheme from ```for step_n in 1:steps``` to ```for step_n in eachindex(step_array) where step_array = [1:steps;]``` reduced allocations by 1. Further changing this scheme to a zeros-array that is steps long reduced allocations by 51. Testing shows that allocations scale linearly with the number of steps, which is incorrect.

So my old method will allocate several intermediate array slices on the heap. If I run away from the inner for loop for each particle, I drop allocations by half with no further optimization. Each operation allocates a new intermediate array. So we have to use the staticarray to get over this problem

``` 
 28.200 μs (3057 allocations: 122.11 KiB) with mutable Vector
 13.600 μs (57 allocations: 78.73 KiB) static vector
```
finally. Well not yet, because we still have allocation scaling with step length, so changing datatypes is no sure-fire solution.

After testing, a problem is shown with how Julia processes complex expressions, such as the position-update:
```julia
position[i] += velocity[i] * stepwidth .+ force_currentstep[i] ./ mass[i] .* stepwidth^2/2
```
Several arrays are here allocated and deallocated within the string, and when we have for-each-step and for-each-particle, we end up with absurd allocation figures and runtimes. When we decompose this expression and we allow the compiler to run  on the entire block of particles positions, we obtain extreme speed and constant allocation for simulation size. The only problem is *readability*. The position calculation becomes:
```julia
positionIntermediate1 = velocity
dumloop_multiply!(positionIntermediate1, stepwidth)
positionIntermediate2 = force_currentstep 
dumloop_divide!(positionIntermediate2, mass) 
dumloop_multiply!(positionIntermediate1, stepwidthSqrdHalf)
dumloop_add!(position, positionIntermediate1)  
dumloop_add!(position, positionIntermediate2)

where

function dumloop_multiply(vec::vectorOfVectors, vec2::vectorOfVectors)
    for i in eachindex(vec)
        vec[i] .*= vec2
    end
end
```
We avoid any allocations within for-each-step by writing over intermediate values  at each step. These are generated ahead of the simulation loop. I've made multiple attempts with nested broadcasting or broadcasted mapping, and these methods result in allocations for each particle, but reducing each value transformation to either a direct substitution or a loop'ed broadcast worked best. Here is where things get really confusing. As part of updating position and velocity due to only the velocity and apparent forces, we have to evaulate 9 for-each-particle loops at every simulation step. So what if we reduced the loops into a single loop which defines a block of numbers to crunch for every particle. That was my original idea in this simulation. Julia hates this idea:
```julia
for step_n in eachindex(steps_array)

        for i in eachindex(objectindex)

            positionIntermediate1[i] = velocity[i] 
            positionIntermediate1[i] .*= stepwidth
            positionIntermediate2[i] = force_currentstep[i] 
            positionIntermediate2[i] ./= mass[i] 
            positionIntermediate2[i] .*= stepwidthSqrdHalf
            position[i] .+= positionIntermediate1[i] 
            position[i] .+= positionIntermediate2[i]

            velocityIntermediate1[i] .= force_currentstep[i] .* force_nextstep[i]
            velocityIntermediate1[i] ./= mass[i]
            velocityIntermediate1[i] .*= stepwidthHalf
            velocity[i] .+= velocityIntermediate1[i]

        end

end
```
Where numbers are crunched, zero allocations occur. But through the nested looping, an obscene number of allocations are made. 
```julia
#TODO make this into a table please thanks
#1. which package? I know in my hunting there is a package big on the tables. Might jsut be Documenter.jl or an extension
```
5 steps, 100 particles
9 loops: 24.600 μs (415 allocations: 16.70 KiB)
1 loop: (9915 allocations: 313.50 KiB)

5 steps, 100 particles
9 loops: 111.600 μs (415 allocations: 16.75 KiB)
1 loop: 29.068 ms (95415 allocations: 2.91 MiB)
For reference, the naive position expression is at least 3x slower the 9-loop on run time, but about half the allocations as the 1-loop.

I am rather confused about how Julia optimizes these 2 methods. I prefer the second method because I can parallelize it infinitely with a macro tacked in front of the for-each-particle loop, but as it is presently written, Julia cannnot handle it well. And I am certain there exists a simple routine to parallelize the 9-loop method.


Another point of consideration is using StaticVectors over MutableVectors to contain dimensional data. We can just avoid loop hell, keep the syntax clean-ish, and improve performance and allocations. But static vectors aren't as workable. You cannot iterate, broadcast, map over them, so then I have no idea how to institute boundary conditions if I cannot analyze or affect any particular value. However, I could instead use StaticArrays methods to affect test mutable-vectors, and depending on the values within these test vectors, I can apply action vectors over the elements of our static vectors. This could make for a beautiful evaluation, and allow me to get rid of the naive boundary function:
```julia
function boundary_reflect!(ithCoord, ithVelo, collector::Collector)
    # can this be evaluated more efficiently?
    #restructureing would allow a simple forloop
    if collector.min_xDim > ithCoord[1] 
        ithVelo[1] = -ithVelo[1] 
        ithCoord[1] = collector.min_xDim
    end
    if collector.max_xDim < ithCoord[1] 
        ithVelo[1] = -ithVelo[1] 
        ithCoord[1] = collector.max_xDim
    end

    if collector.min_yDim > ithCoord[2] 
        ithVelo[2] = -ithVelo[2] 
        ithCoord[2] = collector.min_yDim
    end
    if collector.max_yDim < ithCoord[2] 
        ithVelo[2] = -ithVelo[2] 
        ithCoord[2] = collector.max_yDim
    end

    if collector.min_zDim > ithCoord[3] 
        ithVelo[3] = -ithVelo[3] 
        ithCoord[3] = collector.min_zDim
    end
    if collector.max_zDim < ithCoord[3] 
        ithVelo[3] = -ithVelo[3] 
        ithCoord[3] = collector.max_zDim
    end
end
```
I might have to get rid of this naive function, as it adds 1 allocation per step, but not internally, only in the context of the for each step loop.

At the moment, my attempt to make a ```simulate_unified!()``` for the SVector hasn't worked, it allocates scaling with sim duration and becomes increasingly slow with duration. Some more time spent looking at the work may help, such as trying to restream the functions into a single task, rather than the current split method we are rocking with.

I spent some time trying to understand how we can reduce the 6 loops that calculate position by using broadcasting to achieve *syntactic loop fusion*, but the compiler has no idea how to achieve this in my context and it just continuously allocates temporary arrays instead of looping around each piece. And it is shown that this effort at loop reduction simply recreates the problem in ```simulate_oneloop!()```, which I presently suspect is some violation of cache locality. Here is the undeveloped code:
```julia
        for i in eachindex(positionIntermediate1)
            #positionIntermediate1[i] .= stepwidth .* velocity[i]
            positionIntermediate1[i] = broadcast(x -> x * stepwidth, velocity[i])
        end
        for i in eachindex(positionIntermediate2)

            positionIntermediate2[i] = broadcast(x -> x / mass[i] * stepwidthSqrdHalf,  force_currentstep[i])
        end

        for i in eachindex(positionIntermediate1)
            positionIntermediate1[i] .+= positionIntermediate2[i]
            
            #broadcast!(x -> x + positionIntermediate1[i] + positionIntermediate2[i], position[i], position[i])
        end
        #for i in eachindex(position)
            #position[i] .+= positionIntermediate1[i]
        #end
```

I'm not certain we can read each value of ```force_currentstep```, multiply it by a constant ```stepwidth```, and write the result over the corresponding value of ```positionIntermediate``` without allocation, but we certainly cannot perform anything more complicated while preventing the compiler from creating temporary arrays. So until a better method is found or formulated for this data structure, we may very well have to call a ```dumloop!()``` to calculate each and every modification of our ```Vec3D```, which is going to get very annoying, unless I can embed the function into an operator sign.

## 5. 11 June - Starting up again
My progress has halted for a few weeks as I have followed the ADHD fixation of hyperoptimizing just the Verlet Velocity scheme until failure. It's mostly been really enjoyable, and informative, but it causes problems! Static vectors are far more performant and rewrite prone than mutable vectors, which are allocation prone. But their lack of indexability makes my Naive algorithms more complicated. Naive boundary condition? Forces? Velocity rescaling for energy-maintenance? Dihedral interactions? Virtually anything else? All can use a nested for-loop for each particle for each directional component. So while I have a single-threaded, CPU based simulation engine that outperforms Molly.jl for calculating position and velocity from the velocity and forces that act on them, I don't have much else. I have a force function that must be rewritten, I have a neighbor finder that is innappropriate for this generalized simulation use case, and I have a bunch of slop that needs washing from my source. All of this to say, I have sinned in writing software, and this has been my admission. So now I want to take the opportunity to work on the organization of these docs to reflect what I have worked on and what direction is next. 


But continuing the wonders of optimization anyway: Taking on the method of Molly.jl for position and velocity verlet calculations and adding a generator expression makes it faster than my looped hell expression, bu much more data intensive. Additionally, using mutable vectors with the Molly generator expression is actually faster and more memory volume efficient than the static vector. Now right properly <i>nothing</i> makes sense. However for memory AND allocation efficiency, only the looped-hell expression is winning. (Turns out I just needed some organization and reorientation, with that hyperfixation wobbling.) I still have no idea how to escape this new scenario of either speed or memory efficiency, so we will split simulate along the different methods and continue development from there, allowing for strongly Naive implementations with MVectors, andd maybe eventually we will figure things out to make the SVectors preferred. But I mostly want to stick with MVectors until there is a presented need for SVectors' speed, which currently? Not even close. But maybe once we have some forces, then we will be less confident!