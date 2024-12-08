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
```
_nc = floor.(Int, (xmax .- xmin) / (cutoff / lcell)) 
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

## 29 June 
I extended the dumloop methods to cover the cases for static and mutable vectors and could get a very simple test bank running. The generator expressions for Mutable vecs are worse than loop hell, 7x slower and 38x more data. For SVecs, generator expressions are 30% slower, while keeping about the same data-efficiency as MVecs. Curiously, SVecs take more data the MVecs, but in both algorithms accumulate the same number of allocations. For simple speed purposes, the SVec generator expression is more than acceptable, being 2x faster than MVec loop hell, even with 24 thousand *more* allocations. The ultimate speed appears to be loop exclusion SVec, but these, as previously covered, still bring a tantrum for usability. I hope it improves, but I have no idea how. But I think the compromise could be any indexing functions for the SVec data type just copy! to an MVec and then back to SVec after functioning. If that works! And it does, but would require even more fluff outside of the time-step loop. So mutable vectors for now, since we are losing less than 3x performance, and the copy! operations would probably capture some of that improvement from cache-discontinuities. Maybe in the future I will discover how to make this work without


The boundary reflect. After gutting the entire interior of the function, it still allocates once per step. After type annotating every neatly named variable, still. After type annotating the type of collector to be used in the function call, still! Now it is time to test with other functions, to figure out where this alloc comes from. Getting rid of the function overlay and just tossing the code of boundary_reflect! into the body seems to kill performance by a profound lot. *Some minutes later.*A bug? I created a function no_allocs!(position, velocity, collector) with type annotations that did nothing. I placed this below a commented boundary_reflect!() in simulate!(), and I get an expected number of allocations. But then, I place the analysis code into no_alloc!(), run the code, and now it is allocation hell. Then I remove the analysis code completely, and we are still in allocation hell. Restarting Julia does end allocation hell when no_alloc!() is empty prior to NaiveDynamics precompilation. This certainly makes testing the problem more difficult, but I noticed a doulbe defintion of the same function same method somewhere else, so perhaps fixing that will mean Revise.jl works a little bit better. Curiously, the allocations remain after removing the code block until REPL restart, but execution times drop, suggesting that hey, some new code is no longer being executed (consistent 480 to 415 ms). It seems to me a bug to have evil allocations that appear to be unrelated to the execution of my code. 

Well, regardless of my feelings, the allocation vanishes when the if statement is prevent from evaluating. Provided that any of the if-statements evaluate out to true, then there will be 1 allocation of 96 bytes that is then rewritten across the expression. This seemingly is deallocated at the end of the forloop, and we begin again in the next step of the simulation. However, I am having trouble reproducing this in a fresh file with no dependencies. But regardless of what is being evaluated in context, we still get this allocation. Curiously, if I cleave everything out of boundary_reflect!(), precompile NaiveDynamics, and run@time where the function is called, then the first run of @time generates only 96 bytes of allocated data, but with all the gooey logic, the allocations are much higher. More curiously, replacing @btime causes the allocation to "disappear". Which is to say, the machinery of @btime sitting at the function call is incapable of catching this allocation, when it catches it at the simulate!() function call, while @time has no problem catching the issue wherever it is called.

My next idea stands around syntax, where having to consider either ```position[each][1]``` or ```velocity[each][1]```, that is to say, the float value within the vectors of vectors, results in a 96 bit allocation that is deallocated at the end of the function call. This I now declare follows in the tea-leaf reading of my @profview_allocs that shows an Unknown_type garbage collection allocation occuring in the function call / function space of boundary_reflect!(). Holy banana. My type declarations crusade earlier was the correct direction, but I failed to specify the type of the inner-value. IT IS SOLVED. Now the problem for the code context, how to grab this T. In the user input, we select as a field of the Collector struct for our precision, but just setting T = collector.T inside  simualte!() does not clarify the type for boundary_reflect!(). And now with our types fully allocated, we are down to 19 total allocations, rather than 419. What a good success I think, but this leads us to the next question. I want to have this type stability avoided entirely, but having the position velocity etc. defintions thrown in simulate!() just be renamed entities from our datapaths. How do we get type stability by creating shorthand names or *pointers* to the specific data we are asking for? I am not certain, Molly gets around it by just having less nested structures or by pulling values into the simulate!() function via another function. So we have cured our type instability, but also broke simulate!() to anything besides a vector of mutable vectors of float64s.

Exploring my implementation of the boundary_reflect!(), I realize that a particle that leaves the box and has its component velocity sign flipped may become trapped outside the box if it cannot re-enter the box before the next boundary_reflect!() function call. Instead, the wall should be a secondary force like Lennard-Jones that opposes pasrticle traversal through the wall. The question then is, how do I write a force that will, in effect, flip the sign, but prevent the creation of unrealistic forces? THe simplest solution would seem to be pulling the position back to a small value ahead of the minimum value, which I originally had.


Managing these package extensions so that I can avoid precopmiling Makie until and unless i want to use a Makie feature is very frustrating. [Here](https://discourse.julialang.org/t/should-we-define-new-functions-structs-in-an-extension/103361) a discussion focuses package extensions on only extending functions defined in the base of the package. So I suppose I could provide a dummy record_video() in Simulator.jl, but this wouldn't be extension. This would be substitution functionality. And in [this pull request](https://github.com/JuliaLang/Pkg.jl/issues/3641), the issue seems to be ongoing where people want to modify the methods of 'optional' packages. So manual insertion it is, then.


On debugging NaiveMakie
 we have iteration working correctly, as we can print positionsToPlot for each step of
     positions::Vector{Vector{MVector{3, Float32}}} / each step of recorded simulation
 but these new values are not being plotted
 we changed simulate!() so that instead of modifying copies of the fields of ObjectCollection
     we are moddying the struct directly
 the positions are the same between frames
 new position values are not emerging, in either the loop fusion or dumloop hell method
 I am hoping that changing the collection structure to mutable will magically fix the issue, but i also worry there are some hyjinx related to variable scope, the same stuff taht prevented the chunk indexing from working.

 Delaying the progress of work substantially is the lacking of a dedicated testing package where tests are just in function in line print statemetns that clutter up the place. And as I continue develoing and rewriting, new breakages occur, such as when I cleaned up the stack of equal signs in simulate_bravado!(). Also I have this nasty collection of methods, simulate, simulate bravado, fused bravdo, mvec, svec, and oneloop. I achieved cleaning up the logger mess, but nonetheless it is difficult for working efficiently and easily to go around this boneyard of cold code, even with the use of a scratchpad space for testing specific features, a dev directions page, and a diary to help track my thinking in real time. It is very dangerous to get distracted and forget that I am modifying code in simulate_bravado!() instead of simulate!()

 Anyhow, we have made it to the point where the sys.velocity mutates its values correctly, and even in the simChunk, but when the chunk is passed to the simLog, something goes wrong. In a 20 step simulation, we test the following code inside record_video:
```julia
    println(first(simLog).currentstep)
    println(first(simLog).position)
    println(first(simLog).position==last(simLog).position)
    println(last(simLog).position)
    println(last(simLog).currentstep)
    println(simLog[1].position == simLog[2].position)
```
And it shows that each index of simLog, which are objectCollections, are the same, and they are the same regardless of simualtion length. What this shows to me is that we are not correctly modifying the values of our simChunk and are instead just pointing to the value in the objectcollection instead. But similar testing in the chunk itself shows that its values at each index are different, indicating that within a simChunk we have time progression. But then when this chunk is passed into the simLog, each index of simLog is pointing to the same simChunk, so we can atmost get information only at the last moment of the simulation, it would appear. And this makes sense, as sometimes I see some particles are missing from the mp4. What is strange is that simChunk is just a fixed width vector of objectCollections, whereas simLog is just a multiple of simChunk, for however many multiples of chunk_length fit into the simulation duration.  So my confusion is, why are the 10 chunks of simChunk different from each other, when the same chunks in SimLog are different? Well I just retested it, and they are identical! The code, set right before write_chunk!(simLog, simChunk)

```julia
println(simChunk[1].position==simChunk[9].position)
```
evaluates to true, even though I have proven more recently than the variation in my chunk indices the varation in objectcollection.position for each step. This is what I need a dedicated testing module for, to be able to track the results from a given change.  Anyway, the coding problem does make a lot of sense:
```simChunk[chunk_index] = objectCollection```
is just stating that simChunk at the chunk_index is a pointer to the object collection.
Although according to the docs, performing a[index] = b  will mutate a. So where forces are updated and positionIntermediate1 is upated at each step, it may do well enough to have a for loop or nested for loop to get values updating. So when I print sys.position at the end of each sim step, the values are changing. But each chunk block is not changing. Maybe this is the trouble I read about somehwere of having arrays of structs instead of structs of arrays, Array[Struct] doesnt work the same as Array[Number]. I would most prefer to use ```fill!(simChunk, objectcollection)```, however there is no method matching for this. I have to extend the fill! function to my type, I wonder how much effort that would require. . . Very likely a lot. As fill!() is intended to  be used on a whole array, or part of an array, with a numeric value. The function already collapses in the case of ```fill!(Vector{MVector{3,Float}}[1], MVector{3,Float})```. I require the object I am using to fill an array to be a number, not an entire struct. This is an unfortunate position wherein the only method that appears to work well is the one I have spent a lot of time trying to get rid of, where to store the entire object collection I was pushing each individual field to a parent structure that had nearly the same shape, except each field as a vector of the fields of objectCollection. That method was miserably slow! But I may be able to use chunking with ``` fill!(), as a structure of arrays of arrays will bring me closer to the data points than an array of structures of arrays, and chunking (and more features) should help a lot to speed up time. But for now, the new method is to push the current objectcollection into the simLog at each step. Hopefully I can plan out more how to implement it correctly and be certain to write my notes as i develop.

Back to the video recording, now I have my positions getting printed out directly and I am learning 2 things. The velocities are trending towards zero arbtirarily. No idea why. And the mp4 is only  recording the last frame, as confirmed by visual inspection of the dots and the printed corrdinates. Though, over repeated tries I am not confident in visual inspection, or my changes to the code have caused a switch to the first frame rather than the last frame.

I do not know why GLMakie isnt yet rendering the right video. It is quite the delightful tool, and I am likely less than 40 characters from making it right, but it's function seems consistent with Molly.jl's implementation, and it doesnt seem inconsistent with teh documentation, I just have to provide a number of frames and that at each frame, we update our observable with new values. Maybe i could try a for loop without observable? I think it fell apart at the onset then. Can't remember, haha!!

The record function from Makie is correctly iterating, as an index whose value is printed and incremented inside the function does display correctly. But all the positions are the same, as confirmed inside record, right after simulate. They vary at the end of each sim step, and the values in the log are the same values  the last value of position in the simulation. Something is going wrong with the pushing of the object collection into the simlog, but I do not know how to diagnose it. So I am going to instead push the position to a position log at each step, adjust the return, and adjust NaiveMakie. I have made the adjustments, so now simulate!() returns poslog, a vector of positions where each element is the position at step_n, and these values are all the same! Looking closer, eash time we push!() sys.position into poslog, previously pushed values are updated to the value of the latest addition. That is crazy! I theorize that this is the result of push! being array values oriented, not array of arrays arrays oriented. Append does not work, what about vcat? No, we have to instead push copy.(sys.position). That's truly mad!!! AHHH it works!!!!!!

## July 4, Methods towards improving the development process
1. Develop a testing module to hit the code correctness so that I can change how something works, and know almost automatically whether my changes mean that my project no longer works. Instead of having to develop and test in line if it works right,  then remove the tests and rely on memory that it went okay, and reintegrate the tests when I refactor a given function.
2. Outline the function and motive of a given feature add-in refactor in the dev-directions, then come to the diary and write out in more detail what I want to do, how it's going, what is causing a particular implementation to fail, and what I had to ultimately do to get it right.
3. Have a handy helper sheet that is language specific for important points.
4. Make commits more often, I have been sitting on changes of changes of changes for a while at the present commit.
5. there were others but I forget.



Separately, it is time to put the Velocity Verlet optimization to rest, once and for all. I began to suspect that my loop hell testing was only fast because it helped to get around type insecurity. So I tested:
```julia
# loop hell
positionIntermediate1 = copy.(sys.velocity)
dumloop_multiply!(positionIntermediate1, spec.stepwidth)
positionIntermediate2 = copy.(sys.force) 
dumloop_multiply!(positionIntermediate2, inverse_mass)
dumloop_multiply!(positionIntermediate1, stepwidthSqrdHalf)
dumloop_add!(sys.position, positionIntermediate1)
dumloop_add!(sys.position, positionIntermediate2)

# simple expressions
accels_t = sys.force ./ sys.mass 
sys.position += sys.velocity .* spec.stepwidth .+ ((accels_t .* spec.stepwidth ^ 2) ./ 2)

# generator expressions
(accels_t[i] = sys.force[i] ./ sys.mass[i] for i in eachindex(sys.force))
(sys.position[i] += sys.velocity[i] .* spec.stepwidth[i] .+ ((accels_t[i] .* spec.stepwidth ^ 2) ./ 2) for i in eachindex(sys.position))

```

And the results were strong: at 40 particles and 4000 steps without a neighbor list calculation (so they should all be doing the same number of calculations)
  178.058 ms (8319170 allocations: 317.56 MiB) = dumloop
  541.077 ms (20239170 allocations: 772.27 MiB) = simple expressions
  59.114 ms (5039171 allocations: 191.83 MiB) = generator expressions
I am pleased with this result, but wow the generator expressions are quicker, and at this juncture, I do not know why.
I understand it being faster, but also more memory efficient? I am suspicious. ~~It's possible it is to do with temporary allocating in the copy.() calls that occur with dumloop intermediates that the generators don't have to explicitly state(that I think can be optimized away with a for loop and a fill!()), but consider: The generator expression adds a single in-line for loop, the loop hell is far far more complicated and annoying. Let the compiler do it for me.~~ The generator expressions calculate accels_t and dt but because they are introduced in the for loop they are local to its scope and nothing happens. Fixing that, I still am not seeing the position values change, and I don't understand why. The generator expressions are not working anymore, when they were working earlier without modification. Strange-world we live in. Maybe I had my function calls confused or something of the sort, I swear I saw a video from the generator expression simulate!() that worked correctly. Well, tossing out the the generator expressions that did less than they should have, now loop hell is once again in the lead for superior allocation performance, a 20x perf advantage and a <3x overall memory advantage. And I confirm both methods give a valid mp4.

Moving forward, we can try reimplementing my older scheme to reduce temporary allocations, the one that maximized loop fusion but ended up with performance as bad as the one loop method. And we can also update one loop as well. And while reimplementing maximal loop fusion, I found that miscellaneous allocations were being caused by the ```GenericSpec```s, a type conflict in which the specs were offered only the abstract integer form, and preventing the automagic syntactic loop fusion I ahve been chasing this entire time. I run zero allocations with the following expression, now that *every* variable has a concrete type at or before compile time:
```julia 
        for i in eachindex(sys.position)
            sys.position[i] .+= sys.velocity[i] .* spec.stepwidth .+ ((accels_t[i] .* spec.stepwidth ^ 2) ./ 2)
        end
```
The compiler has shown that it can get around the need for explicity intermediates, reducing data generation and at least matching the speed of the dumloop method. And it should be simple to exceed the speed with this fully compressed expression.


We are working towards preventing particles from escaping the box. We readded position locking, which only works if the position coordinate doesnt NaN, as NaN is not greater or less tahan any number. we can increase the duration of the simulation and decrease the velocity range, to hopefully increase the probability that two particles do not appear extremely close together and then fly off at incredible or overwhelming speeds. Velocity rescaling could also help catch some particles, as we aim towards conserving energy in the system. I think tuning down the sigma and epsilon terms could also help, whereas before tuning them but not initial velocity was not enough?

Also, at 4 million steps and all the nonsense in my background, we are hitting a major memory wall to the point this program may crash.

## July 18, trying to fix NaN and Inf
The forces are not being adequately zeroed out, testing shows that no matter how I apply the zero() function to my datatype, nothing indeed happens to the interior values of the argument. But this works in a generically prepared array of MVecs on my scratchpad. A for loop and map routine is able to eventually solve the problem, but it's syntactically ridiculous to have:
```julia
    for each in eachindex(force)
        map!(x->x, force[each], MVector{3, T}(0.0, 0.0, 0.0))
    end
```
instead of
```julia
zero(force)
```
So a new task for the devdirections, extending the zero function a Vec3D? Though maybe the zero function is just a wrapper for a fill with zeros function.


Towards velocity rescaling, life appeared to be working fine until the rescale_velocity!() function found its way. When calculating Ti in the nested loop, Ti turns into NaN at the first calculation. Why is that? 
```julia
Ti += (2/(3 * objectcount * kb)) * v * mass[each]/2
```
When Ti is calculated in this expression for each particle, the presence of v---

hold hold HOLD the phone. whenever two particles get close enough for an Epic amount of repulsion, they go NaN and this somehow breaks every other particle

## July 31, a fun case study
```julia
function threshold_pairs(list, threshold::T) where T
    
    return [list[i] for i in eachindex(list) if list[i][6] ≤ threshold]

end

function threshold_pairs_old(list, threshold::T) where T
    thresh_list::Vector{Tuple{Int64, Int64, T, T, T, T}} = [] # this is a silly fix
    # would it more perf-efficient to define a threshold list as long as the unique pairs list
    # at small n particles, and just reorder the threshlist between valid and invalid values
    # and jsut instruct functions to use the 'valid' region of the array?
    thresh_list = []
    for i in eachindex(list)
        # replace with named tuple?
        if list[i][6] ≤ threshold
            push!(thresh_list, list[i])
        end
    end

    return thresh_list

end
```
in this case, since the type of thresh_list inside threshold_pairs_old is already secured, the generator expression does not seem to improve speed, and saves up to 2 allocations. But if the type of the array we are push!-ing to inside a loop is insecure, then we get heavy handed allocations. Fixing this insecurity in the unique_pairs function by using a generator expression and not with type-annotations offered an impressive speed up: 912.800 μs (19998 allocations: 481.22 KiB) -> 17.200 μs (7 allocations: 86.34 KiB), which greatly speeds up the initialization step.


### NaN and Inf are unavoidable
and statistical consequences of having extremely powerful forces involved for a given simulation temperature, box size, and temporal resolution. We can only reduce the probability of occurrence, or maybbe drastically increasing temporal resolution so  that a LJ force can balance repulsion with coulombic and LJ attractions, rather than generating some obscene repulsion.


## 7 August, bitmadness
```julia
# this is specifically disallowed in Julia
a = "000000000000000000000000000000000000000000000000000000000000000"
println(a[1])
```


from stackoverflow to set an nth' bit to 1, 
```julia
N |= 1<<(n-1)
#where N is a number and n is the nth bit we are flipping
# or
#but 2^n is lsower than 1<<n
```

Where
```julia
N ⊻= 1<<(n-1)
``` 
will flip the state of the bit at 'index' n.

## 10 August, morton coding towards a ray castable BVH
Here is to spelling out the process as I am getting lost in the sauce.

Given an array of positions, we must
1. Define boundaries of an array of AABBs
2. Split the simulation window into a user selectable set of bins
3. Sort each dimension least to greatest
For the morton interleacing, interalce from the centroid or bdounaries?
4. Assign the new ordering an index ranking for each dimension
5. Interleave in a bitwise fashion the first 10 bits of each index ranking into a single Morton code. I don't know if there should be 1, 2, or 3 morton codes for the centroid, the min boundary, and the max boundary. I guess we will find out later in traversal!
Or we find out now, the morton codes would alll be identical bc we end up just having a grid based around *these* values or *those* values which are related to each other by a constant and universal shift value of the critical radius.


Oh this is annoying, given the constructor:
```julia
bvhspec = SpheresBVHSpecs(; floattype=Float32, 
                            interaction_distance=0.1, 
                            atoms_count=length(myCollection.position), 
                            bins_count=length(myCollection.position) )
build_bvh(myCollection.position, bvhspec, myCollector )
```
I am getting a no method errror

```julia
ERROR: MethodError: no method matching build_bvh(::Vector{StaticArraysCore.MVector{3, Float32}}, ::SpheresBVHSpecs{Float32}, ::GenericRandomCollector{Float32})

Closest candidates are:
  build_bvh(::Array{StaticArraysCore.MVector{3, T}, 1}, ::NaiveDynamics.SpheresBVHSpecs, ::Collector) where T

```
which is very annoying, what part of this have I done wrong? Curiously, the REPL had highlighted ```::NaiveDynamics.SpheresBVHSpecs``` red, causing me to think that was the problem. Upon the reload, it became clear the problem was to do with the Collector typing, as subtype where a field-less abstract type was expected. A similar sort of error appears later on in the track. So if we close Julia and load it again, is life now fixed? Yep! As far as I can tell I changed nothing, but reloaded the REPL recursively to ilet ti fix its own. This language TRULY is not meant to be used with a text editor other than its own command line.



Next!
```julia
Failed to precompile NaiveDynamics [ef6c0610-87db-4406-b8dc-01afecc28d92] to "C:\\Users\\kucer\\.julia\\compiled\\v1.10\\NaiveDynamics\\jl_FB4C.tmp".
ERROR: LoadError: UndefVarError: `update_bvh!` not defined
```
Truly confusing! I have no idea why this doesn't work. Curiously, this error changed from covering my build_bvh function to my update_bvh! function just when I changed their order. Ah, this one was my fault, I forgot a comma in the previous line.

Onwards!
```julia
struct GridKey{T} <: AABBGridKey
    index::T
    morton_code::T
end

ERROR: setfield!: immutable struct of type GridKey cannot be changed
```
If we are running bitwise arithmetic on the morton_code field of an array of GridKeys, should this error out? I don't think so as we are simply manipulating the bits of the value, not changing its type nor changing the shape of the structure. But who knows! Mutable!


Sometimes when evaluating the length of the common prefix between L of i with L of i-1 or L of i+1, both turn out to have a common prefix length of zero. How does this happen and what does it  indicate? It perhaps indicates that the internal node, I, has only 1 child leaf, L. The question then is, how is this possible and why does it happen virtually every single time? There are other artifacts too in bvh solver at the moment, such as j sometimes being wrongly posed as some value way outside of the range, as currently there is no control circuitry when we finally get a common prefix shorter than the root node to pull l or lmax (wwhichever it is) back to the previous value and work it's way incrementally to find the true maximum index away from i. So in the ordering of our morton codes, I am now curious. It turns out, our ordering is at fault, from least to greatest rather than lexicographical, or left to right, bitwise.  Consider the bit strings of 0, 1, 2. They are ...000, ...001, ...010. Finding which direction "d" should go from the perspective of the morton code = Int(2) is impossible, and it likewise impossible for others. We need to work on our morton sorter! Though this now creats a new problem that inspecting either direction may yield a common prefix greater than from i to j. So now we just have to compare the larger of the two, I suppose?

Also, my methodd has no protections that each morton code is unique. ALSO, this adjustment to the sorting does not resolve the problem, I still see a zero prefix surrounded by 1 prefixes and vice versa. How perilous! Upon further investigation, I now understand the meaning of prefix in this situation. Hah! Fixing the prefix finder and tossing in some control logic for when j = i +/- 1 solves that issue!


On the sorting method, we are in an annoying spot!
```
00000000000000000000000000000000
00000000000000000000000000000000
00000000000000000000000000001101
00000000000000000000000000000011
00000000000000000000000000000011
00000000000000000000000000000110
00000000000000000000000000000111
00000000000000000000000000000111
00000000000000000000000000001000
00000000000000000000000000001000
```
I would assert that this is incorrectly sorted, but maybe bitstring() instead of string() sortby will help. I dont remember why sortby bitstring didn't work, though. Ah, I may have changed from bitstring because sometimes the morton codes are identical either side of the comparison, but that is no fault of the sorting process but a fault of the morton assignment process. There are 2 tricky details. The first and simplest, two particles can exist arbitrarily closer together in real space. The second arise from Karras implying that morton codes are not bitwise interleaved integers, but float values themselves that are interleaved. But from then on in the text and elsewhere they are treated as integers. This difficulty came to me first wwhen I tried implementing the digit interleaving and found I really could not effectively represent the positions of a Float32 as an Int32 but just multipliying the Float by some large factor, and then rounding it and converting to an integer. The problem I ran into was overflow from trying to catch too many digits of the float in the int. 


```MethodError: Cannot `convert` an object of type Tuple{Int64, Int64, GridKey{Int32}, GridKey{Int32}} to an object of type GridKey{Int32}```
So here the problem is in trying to set the Julian version of a pointer, a Ref. But the Ref expects a fixed type, but an internal node may have at the left a leaf, and thus a gridkey type, or an internal node tuple. How difficult! We can get rid of the problem by creating an INode struct:
```julia
mutable struct INode{T} <: NaiveNode
    leaf_indices::Tuple{T, T}
    left::Ref{Union{GridKey, NaiveNode}}
    right::Ref{Union{GridKey, NaiveNode}}
end
```
but this introduces for the first time in this code base the difficult side effect of changing the types of the fields of a struct *at runtime*. This is likely an unavoidable structure of a garbage collected language, which allows me to abstract memory management away, but prevents me from having an array of pointers that purely point to fixed addresses in memory. I suppose in a HPC  C-code that could be terrifyingly bad practice, but how else?


## 24 Sept, I remain silly
1. I want to reconstruct this from a simple stream of thoughts, mostly of 'this is the bizarre problem i have, everything must be borken and i am the victim' to audience oriented blog posts and shift all the stuff for me that goes here into a do not include section of the git(?)
2. i want, maybe separately, a state of the union to sign off from a coding session discussing what i accomplished, maybe some metrics around it, and immediately directly going forawwrd, so as to reduce the time of rediscovery. i leave for a month and come back with no idea if my algorithm generates a hierarchy. I believe what ti does is correct, but it doesnt do anything to make it traversable? as in, how do we determine the traversal tree. or at least, what the hell is Karras using atomic counters for??
3. I also want more formalized decomposition of an article series for my work. an intermediary between the initial reading, the paper decomposition and the final coded product. The intermediary should be where i concretely describe whwat i am trying to do, what my current implementation looks like, questions about the meaning of the text(s) and their interactions. this would better ground the theory, but also is maybe not super practical, as i prefer to learn by running my face into a wall. but maybe it is soomethign to try in order to implement the super method from prokopenko, especilaly as their paper is centered around the exact kind of decomposition i am now thinking of, albeit an abridged version. I think it is important towwards the problems of uninformed and meaningless meandering, misunderstanding the theories and ideas altogether, and misreadings/underreadings of the source text hiding information from myself.
4. i want bettter structuring of the neighbor seearch file, yeesh i am having a hard time ith the redundant data structures, trying to figure out where i should go! i dont want to change to many thigns and bash my head for subtley  breaking many of them.


Starting from ```Inode[1]```, how do I determine which Inodes are the left and right children? We determine boundaries, where ```INode[1]``` is always the left child of the root and ```INode[n]``` the right, and they have bounds as half of the box. So need to add to the data structs our boundary tuples! Two coordinates recombined will grant the 8 points of the box that encloses the 'van der waals' sphere surrounding an atom.

WE construct all INodes in 'parallel', how do we affix a 'parent pointer' to the leaf node at the end of each parallel construction? Running through the logic, if the minim value between the Inode being considered, i, and the Inode beign determined, j, is the same as the breakpoint between them, then the left child of Inode i is the leaf node gamma. (there is the equivalent max case too) Therefore, it follows that the parent of the leaf node gamma is Inode i. or, 
```
if min(i, j) == γ
    left = L[γ]
    L[γ].parent_INode = i
```


So if we had an atomic add for each node, how do we determine the boundaries of each inode in a parallel way? IN the present implementation, I believe I require another for loop after bvh_solver over each leaf node.
1. If right or left child, leaf node will give parent the min or max value
2. Do we ahve to march up the tree to determine the boundaries of the nonterminating Inodes, or can we nab the boundary in parallel? Feels pretty sequential to me


### hopefulyl i will go back and read my work
but this is still broken in bvh_solver in some way or another, as the vsits counter does nto exceed 2 even for 500 particles, and all visit counts past about inode 15 are zero. Many things may still be breoken, new logic and old!


## 25 Sept

I dont understand the sense of Howard 2016's use of do block. It feels like a loop that evaluates forever and it depends on v eventualyl being zero, triggering a return. But we always increment v by 1 before evalulating if it is zero 
I also don't understand we cancel out the first leaf thread, and only allowing the second to evalutate the parent. Double counting feels inevitable. I don't even understand the ssignificance of the visit count. I can appreciate that the least visited nodes are terminating nodes, and the most visited are the root. But we already have the reference unions telling us which data point to go if our particle fits in this or that half of an aabb, whic by the way, what about this halfway marking value???? How do we decide to go to the left or right neighbor if we don't know the barrier plane in space between the left and right children of a parent? Maybe we have to hope julia can handle this internally?? ugh.

And then the merge part, first we decide which boundary is larger or samller between 2 things, then assign

Gods this is messy stuff. I am thining through a 4 parter logic whwere we have to ask what data structure we have and which data structure has the smaller value. That is terrible and I refuse to think this is the awy it happens. But INodes have left or right, and while we can refer to thjose


Okay so, Howard2016's do block is either incorrectly notated or it is not a loop construct




## 26 Sept
so we have parent pointers. How do we have a for loop starting from every leaf travel upwards and give every INode their boundaries? The lit describes some elimination process so that only 1 thread processes 1 node. I have no idea how to achieve that with atomics, and probably don't understand atomics either. But I believe we can make do by reducing our work efficiency. The expression,
```julia
I[p].min = (I[p].left.min < I[p].right.min) * I[p].left.min + (I[p].left.min > I[p].right.min) * I[p].right.min
```
will grant some Inodes with incorrect bounds, for all parents of a leaf and a branch. But the leaf nodes at the bottom of the tree should ideally be processed after these elevated branches are processed, thus hopefully resolving the issue. But this is not necessarily true, and it makes this neighbor finder nondeterministic, which is the point,, moreover, solving it sequentially is bound to create artifacts. For now I suppose we have to expand the work further to catch cases of nonanswers. 

Wait, this expression won't work anyways because the 'min' is a mutable vector, not a single value. Darn! WSwe have to calculate distances anyway to follow this scheme! Useless in place distances just to determine which point is overall closer to the origin than the other. We could adjust and just right the overall sum without squaring by assuming a unit square simulation, like we were supposed to anyways. Hah!


27 Sept
Another query, whty is the Inode procedure not performed on the root node? It is incorrect! But i leave the question here that i might find it later.


The new problem I am faced with is trying to process the root correctly in the BVH solver, so that wew can work on traversal. The root node coveres wwhat range of the simulation, and what are it's children? Previously, I had this set to be itself and the last internal node, which is impossible. Ah, it is not set to itself and the last internal node, but to itself and the last leaf. Hence the name, *leaf* indices. Whoops. The root necessarily covers the entire range, so it's calculation perhaps can be simplified to only work on finding the split value to determine which internal nodes are its children. It is important to keep in mind that the array L is sorted according to its morton codings, so covering the netire range of leaf nodes *means* covering the entire simulation window. However we still have to evaluate the conditionals at the end because a sufficiently small tree can result in a root with a leaf child.


The second problem surrounds  increment and decrement. If INode number 2 is selected, can it only consider backward to INode 1 and not further back starting from the end of the INode array?


When we consider INode 1, wwhat if d=-1? In Karras, there does not appear to be allowance for an exception in the case of i=1 d=-1 (or in the paper's case i=0), nor is index wrapping apparent.

WAIT, what is the consequence of their comment, delta(i,j)=-1 when j is not a member of Inode array?




The struggle now is that Howward 2019 uses stackless, where if a particle overlaps with a node at all, then you proceed to the 'left child', and if you dont you proceed to a rope to the next child. Marrying this method would involve a lot of rewriting of my algorithm outside of the traverse function. So oour code can't look the same! Perhaps we can get on wwwith this stackless train later. For now, let's just get the stackfull approach moving.



And towards the stackful implementation, somehow my hile loop is noww failing to correctly loop and work an INode down into a gridkey. The while loop should only exit after appropriately testing a leaf against an atom or ending up with a value that is a leaf/gridkey. If I am correct, the current set up ending up with a nothing result happens whenever one of the children of a node is a leaf, as there is no mutation of 'n' at that point. I think we are supposed to consider that the termiantion of that 'thread' and move on?


Now wew are up to the next problem in which the neighbor list is written, and then unwritten. And that was much simpler, it was because I wasn't returning 'neighbors' at the end.



Hrm, as constructed I am not confident binary search is being performed here. Suppose 1 sphere overlaps both the left and right children of a node? We condityionally must evaluate both halves. But then doing that will boot loop us into the hell of nonhalting. And e are taking a result and overriding it before doing anything useful with it. This structure is not working, in other words. We would need to do a thread splitting activity to independently solve cases where an atom overlaps ith both children, or otherwise miss some pairs.

So the hypothesis is that stack management here means noticing that hey, this INode has 2 INodes as children and my query overlaps with both of them, now I need to test both of them, if only there were some way I could set up checking book to keep track of all of the places I need to go to next. Even if that is not stack management, it would be complicated and expensive to deal with. So now we look to skip directly to Prokopenko and Lebrun-Grandie 2024 as their implementation is apparently the best but likely the most detailed. Lovely!


Key detail is that the delta-star function changes meaning from computing the common prefix to a 'simpler XOR evaluation'. Apretrei 2014 defines delta(i) as the index of the highest differing bit between the rang of jkeys covered by node i. So let's get testing! If 2 keys are identical, then Prokopenko 'aguments the key with a bit representation of its index'.


Okaay so struggling right at the start. Prokopenko uses del(i), but at the first step of the algorithm when we don't have an assignment for the other index of L to compare a given Morton code to. Apetrei a leaf node to allways have a range of [i, i], but if this is a self reference, if del(i) in literature is to mean del(i, i, L, spec) in my code, it will always and forever be zero. What is the second value I am meant to grab by the thought of del(i)? Found it! Reading the paper and taking notes about ti first is really truly helpful, hopefully I do it sometime! They use delstar, which means del(i, i+1)

And the next problem is adding skip-roping to the code requires that L have a skip value, hich is a reference to what? Upon initialization? God damn it! I have to create it before I refer to it?? If we arrive at a leaf, we always have to know the next place to go. Howard 2019 uses the sign to differentiate between one side and the other. Brilliant. Thank you literature, I am sorry for my insults. There is a sentinel node, which is noted as artificial. So until proven wrong, I shall note it as zero.

AtomicCas is a function that runs atomic compare and swap. It compares the contents of a memory location of its first argument wwith its second argument. If they are the same, then the memory location is overwritten with a new value indicated by the third argument. Thus, ranger = atomic cas(storep, -1, rangel) should mean that if p is -1, then it becomes the value of rangel. The value returned by AtomicCAS is the original in memory value, not the new one. is compared with the initialization value of probably negative 1 to determine if the value was change. 

Next problem! We already have cases of j being out of bounds to dummy value our calculation, are we supposed to do the same for i? I suspect so because we subtract 1 from i for Apetrei and Prokopenko, but I'm not finding strong textual evidence for this. Maybe karras or even Lauterbach have an idea.


Nextr! The Prokopenko scheme requires parallelism inherently, but I am getting a ReadOnlyMemoryError. This was fixed, seemingly, by turning the Apetrei index 'p', into parray for each leaf node. Makes sense, don't want different threads all overridding the same value or otherwise having their own meaningless value.

Nexter! How do we resolve 2 morton codes being identical? In Karras 2012, it is suggested to concatentate the string representations of the morton code and the index of L that a morton code belongs to. But they say 'in practice, there is no need to actually store the augmented keys--it is enough to simply use i and j as a fall back if Li == Lj when evaluating del(i, j).' So instead of calculating that concatenation in a temporary value and hoping the compiler will hold it around for rewrites later, we will just use the index values directly, and hope for goodness.

Nexeter! It turns ourt *none* of the threads are entering the exit condition, meaning the program does not halt. I waws hoping fixing the other things would get me far enough along towards a resolution, but it is not so. The key is that the while loop proceeds in Prokopendo and Lebrun-Grandie until i = 0 or until range-right or range-left receive an unacceptable error. The text describes I0 as the root, which my ordering does not have! Maybe that will help. Checkin up on the left and skip connections, we have a value of -11, which is worrying, there should only be L-1 internal nodes. The great difficulty is that the threads are not in lock-step but instead executing randomly. So let's try to investigate if it will still run without multithreading. OUr exit sequences should work regardless of threading. Threads rarely exit the while loop or hit the return calls, and usually a ranger or a rangel value gets stuck at 1 place and function evaluation does not mutate it.

One fail case arrises when dell == delr, and seemingly another wwhen rangel == p == ranger == q. HOwever, if i == rangel == ranger and p[i] == q, then it seems to exit just fine through the if-return. This is what test driven development is for, i am preparing these tests, albeit slowly, they should go in their own test section so the same code can continue to exist! I do Have a logical inconsistency, Prokopenko and Lebrun-Grandie have an if else nested in an if else, where I have an if elseif else. Let's tighten these up first! A second inconsistency in the same block, del-star (i) means to my code del(i, i+1, L, spec). In the paper, they use del-star(i + 1), which to me should read as del(i+1, i+2), yes? Well, a thread still exits badly if i == rangel == p[i] == q. One question I have a this point is around atomics, I have p as an array for atomics, one for each L, shouldn't rangel and ranger and q also be arrays of atomics to each L? Or what values, if any, should be held as atomic and given unique memory locations for each index of L?

Oh my god, I am missing and entire if statement before the atomic cas. WWHATTTT. No! I am NOT!. I waws checking against algorithm 3. Damnit!!!!!!!!!!! Finding out I was looking against the wrong source is even worse than not correctly reading it the first time. I waws so excited because I thought I had caught my mistake, afterall I was taling about "oh this is a problem" with rangel == ranger. That's frustrating.

Okay, we found a legit typo, one 'ranger' was supposed to be 'rangel'. And! A second instance of it! Even wwith these corrections, it would appear the code is still broken. But I would like to add for the paper's syntax. They have 'repeat until i = 0', which means what? Does it mean evaluate at the word repeat 'is i == 0', or only at the bottom immediately before beginning the loop? If Prokopenko and Lebrun-Grandie intended an evaluation at the top, then they would have just used wwhile. Additionally, this structure means that the first Leaf node will never be processed. So we need to have a while loop which goes on forever, with an escape sequence at it's end. There is an important detail about this, at the level of the for loop, i is the index of the Leaf nodes. But within the repeat until 0 loop, i becomes the index of the INodes. Adding a fix for this presents a different problem, given the stack:
```julia
i = 1
rangel = 1
ranger = 1
dell = -1
delr = 3
p = Base.Threads.Atomic{Int64}(-1)
q = -1
```
we get an error crash out for trying to access I[0]. Rangel was set by the atomiccas to zero. I understand it that this is the first 'thread' (the print statements show this), so the atomic value should be changed, causing the check value to become negative 1, and the thread should exit at the check. Final answer: another typo. When checking if less than 1, we have ranger and rangel, when it should only be the range value that was mutated (or not) in the previous line, (here as rangel). That error is now sorted, back to the looping.

#### Also dont forget that the computation of the bounding boxes is done during the cpmutation of relations!. They just leave that bit as an exercise to the reader


We need to fix morton_coding, we shouldd NOT have redundant codes at 10 atoms. We have a mode of addressing this in del(i,j), but it's still ridiculous for this to not work already.
Where m is the output, and e is the expected code.
```julia
x00000000000000000000000000001010          
y00000000000000000000000000000111
z00000000000000000000000000000001
m00000000000000000000000000001010
e00000000000000000000000001011110

x00000000000000000000000000001000
y00000000000000000000000000000010
z00000000000000000000000000001010
m00000000000000000000000000001010
e00000000000000000000101000110000


x00000000000000000000000000000010
y00000000000000000000000000001001
z00000000000000000000000000000100
m00000000000000000000000000000100
e00000000000000000000010100001010

x00000000000000000000000000000110
y00000000000000000000000000000101
z00000000000000000000000000000110
m00000000000000000000000000000100
e00000000000000000000000111101110
```
I can't pretend fixing these values will save my day, but the actual ansers are *dramatically* different from the generated answers.

So working on a neww fresh method for updating the morton codes. I have spent years down in the mines and I just don't get it. I made a setup to shift bits around:

 if y = 1 m = 1 n = 2
        then << 31
        then >>> 31
    if y = 2 m = 2 n = 5
        then << 30
        then >>> 31
    if y = 3 m = 3 n = 8
        then << 29
        then >>> 31

We ask a number  to left shift 1 pace, it does so by 2, by three, by -- why am I incrementing how much I want to shift the Morton code by. . . . I only want to move it by 1. I spent a while fudging around with the algebra so I wouldn't have to introduce if-statements (and thus pursue a very complicated implementation). But I can run throug the hidden f(x) analytically by iterating over the index of the input and over which dimension I should be picking. With those iterating away, whichever n'th position of the output morton code is doesn't have to be determined, it just has to shift left 1 bit per lowest level iteration. Nice! Now to reverse a bitstring,,,,,,,,,

The future will hopefully directly generate ithout having to generate and then reverse, but currently, just trying to reverse the order in which bits are added to the result changes the result in an as-of-yet unknown way. A function call, reversebit(n::Int32) allegedly modifies a data structure, and it is seen the result that it does not. At first I was going to freak out, but *inspecting the function closer* it is not a mutating function, takes an argument and offers a new result. Ah Hah! Okay so I reversed the outside order, not the inside order also. The sides were good to 6 and 8 bits in, but that's only because those wwewwre symmetric until then.

#### It is very important to brainstorm on a better way to develop code, because writing it out of tree andthen transplanting it into context usually breaks context and requires a lot of annoying syntax adjustment. But the API is not extensible enough to allow me to test each and every function. I think all Big Functions need to be opened up to a public API for in context testing and development.


I am still hunting for logic errors. For a moment I believed I caught one about the use of do - until. It is, I suppose, somewhat ambiguous. Do the authors mean that the entire indented block inside the for loop, includign the repeat block, is to be re-repeated until that block hits internal escape sequences or, at the scope of the entire block, a value is obtained of i = 0 (i=1 in our case)? Or does the 'until' line only bind repetition of the repeat block? Considering the fact that the 'until' line is at the same level of indentation as the repeat, as well as the lines of code preceding the until block, I believe my original coding was correct. How would we show this syntax, in this alternative form? And, as it stands, the code breaks because it becomes stuck inside the while loop. Adding yet another loop around that while loop should have no effect. So we have to look for some other detail, some other piece of evidence. For instance, are the morton codes being sorted lexicographically still? I see that they are being sorted in value order, but what about the bitstring? Changing out the sorting method returns me a sorted morton-coding of
```
1357
188
196
2259
283
362
423
560
```
which certainly feels closer to being lexicographical. Let's try to narroww in what the literature says about sorting. Lauterbach 2008 sorts by increasing value of the Morton codes. Karras 2012 uses unclear language here, stating that they would 'only consider ordered trees, where the children of each node--and consequently  the leaf nodes--are in lexicographical order. This is equivalent to requiring that the sequence of keys is sorted, which enables representing the keys covered by each node as a linear range [i, j]. Using del(i,j) to denote the length of the longest common prefix between keys ki and kj, the ordering implies that del(i`, j`) is greather than or equal to del(i, j) for any i` and j` that are elements of [i, j].'

This can be tested readily, because Karras et al. helpfully lay out the verification procedure for whichever sorting method we should use. However, I am not using Karras's version of del(i,j), I am using Apetrei's. Apetrei does not make it obvious *how* the keys are sorted in their method, neither do Prokopenko and Lebrun-Grandie. Howard 2016 goes into greater detail of their morton sorting routine. They prepend a bitstring for the particle type to the code, then sort by particle type, then sort 'along the Z-order curve'. 


Now we begin to beg wikipedia for more information, as well as analyze the code. My code will either sort the values by increasing value, or the lexicography of their integers, *not* the integer bitstrings. But I think there is a fashion of correcting this. I have some confidence our sorting problem is the largest culprit to progress in this algorithm series. The new challenge is trying to trick sort!() into doing what we ask.

By pulling the reverse bit method back from the Leetcode, it finally produced for us a lexicographically sorted morton code array, but our tests for any interanl del being at least as large as the boundaries del returned false. If it worked, we would have beeen able to return to the original morton encoder, where the order is generated reverseed from iterating upwards instead of downwards.

Just by running and rerunning, we randomly got a series with the desired del() values between morton codes. HOwever, this grouping still resulted in several bad thread exits where they were able to NOw I have to specify because I realize I have been fudghing details here.
```
sort!(L, by = x -> bitstring(x.morton_code))
1.  will sort morton codes of L by left to right, so essentially the code with the most 
    left zeros (that which has the least value) will rise to the top of L
2.  never grants del(i`,j`) >= del(i,j)
3.  some threads randomly exit well

sort!(L, by = x -> bitstring(reverse_bit(x.morton_code)))
1.  will sort by right to left, so most zeros on the right side of the morton code
2.  randomly grants del(i`,j`) >= del(i,j)
3.  some threads randomly exit well

sort!(L, by = x -> string(reverse_bit(x.morton_code)))
1.  sort lexicographically (left to right) by the digits of the reversed integer
2.  randomly grants true evaluations
3.  some threads randomly exit well

sort!(L, by = x -> string((x.morton_code)))
1.  sort lexicographically (left to right) by the digits of the original integer
2.  randomly grants true evaluations
3.  some threads randomly exit well
```
We should not favor results in which true evaluations are guaranteed, if that can be engineered. That sense of the common prefix is limited to Karras, and is algorithmically unnecessary. I believe normal lexicographical sort is the most correct method of morton sorting., so the first one in the above block. So we must look either upstream or downstream of Morton sorting for our problem.

We move on to the next phase of testing.
```
running at least 5 informal tests:
At 1 object, thread exits well
At 2 objects, threads exit well
At 3 objects, 2 threads exit well, 1 does not
At 4 objects, 4 threads exit well
At 5 objects, 3 or 4 threads exit well
At 6 objects, 3 or 4 threads exit well
at 7 objects, 3 or 4 threads exit well
    thread 2 sometimes fails, 4, 5, 6 always fail
at 8 objects, 5 or 6 threads exit well
    there are several 'modes' for where threads fail, unlike all previous results
beyond 8 objects, the ratio of passing to failing threads decreases.
```


Also, it is informally noted that the thread which fails to exit well does not appear to change position, it is held to the index of L
that if it has that particular index, it will not exit. And in perhaps everycase, the first thread exits normally.

Next experiment, somethign about my implementation of the CAS behavior is incorrect. Because if we have a 3 atom system, the second atom/iteration of the for loop will always fail. The first iteration will change the atomic value, become an invalid value, and then the thread will be closed by the following if statement. The second iteration will run the cas again with the changed result.

AtomicCAS compares the contents of the first argument with the second argument. If the first and second arguments are the same, then the first argument is rewritten with the third argument. It returns the value of the first argument prior to rewriting it. If they are different, then it appears that nothing happens. The effect of using this is that the first thread to encounter this code block will be closed out, whereas the second thread will be allowed to proceed on.

But there is something nagging at me. Here is one detail, there should be a store array equal in length to L with each value set to 1. The indices of store are referenced by p. So wwe messed up the logic around the p atomic, not the atomicCAS. Hooray, we have fixed the looping evaluation.

So now the problem is when evaluating i = 1, we get down to asking p to be 0, and we are trying to index the zero'th place of our p-array. Suppose we just try to add one. It is essentially that same burning question about the wrap-around. In Prokopenko and Lebrun-Grandie, following their pseudocode logic, if i = 0 and dell > delr, p will be set to -1. And the algorithm will ask to loook up index -1. Handling this without custom if statements here in Julia will be truly annoying. ArborX, harboring their implementation, is a C++ clang-16 format. And according to a quick Google, negative indexing in C++ does not wrap around to the end of the array but starts looking atmemory locations in front of the array. So, what? what now?? I tried shifting the access by 1, so p is ranger + 1 or rangel, which works because ranger comes from i, and i will be at most 1 - length(L). Now the loop is broken and stacklessbottom_bvh can execute to its heart's delight. Onto the next problem, are all of the left and skip nodes being generated accordingly?

None of the left children of L have been changed, and an assortment of L.left, I.left, and I.skip have not been changed from the invalid value of -length(L). Threading the for loop does not have a great effect. Remembering myself, L.left should always be a null value. I am going to change these values back to zero instead of some impossible negative number, because a serious portion of Internal nodes should skip to the sentinel (marked by 0), wwhen only 1 or 2 leaf nodes should skip to the sentinel. I say 1 or 2 because the last i == length(L) cannot be considered at all, whereas i = n - 1 can be considered, but is always caught by the sentinel evaluation. But I still worry, there are several internal nodes that have both children as sentinels, which is definitely *not* right. I mean, to the credit, looking at 3 skip ropes of L, they dont point to these no nothing I nodes. Many leaves point to each other. Which is perhaps more maddening than anything else. If there are no complaints, now it is time to calculate bounding areas. 

The issue is, if we hack p to work nicely, we don't know enough about the behavior of this code to determine why it is still failing. It is plainly obvious why it is erroring out, if the code of Algorithm 4 from Prop and Leb is implemented, it will crash out immediately by trying to access the store at index 0. After days of working on this and studying every part, this feels like a brick wall to every other problem I have encountered.

Options at winning a functional and surviving method
1. I could try to add a series of evaluations that given a negative number or a value too large, wwrap around to the otherside and idnex backwards from there.
2. Find where my code's problem is in relatino to the text.
3. Find and interpet the ArborX code, which doubtless will be mcuh improved by now. Find it, chronicle it's hiding places, and figure out how they do it and compare back to here.
    looking on tree construction.hpp it looks like the p = ranger is p = ranger + 1, and it is p = rangel - 1
4. ?


#### No wait the complaints are starting up already, traversal starts at I[1], and in every case so far, I[1] is empty.
1. trying morton codes
    sort(string(L)) does not fix t he problem and introduces many many redundant references.
    sort!(L, by = x -> x.morton_code) same as above
    sort!(L, by = x -> bitstring(reverse_bit(x.morton_code))) same
    sort!(L, by = x -> string(reverse_bit(x.morton_code))) same

## 8 Oct.
So after digging through ArborX, it was found that their delta-star returns the typemax for out of array bounds of index values. So now the program will run without issue, leaving us in the wastelands of making the algorithm run correctly:
1. The for loop is only supposed to run across the number of Leaf nodes - 1, meaning that one leaf node is entirely unprocessed.Is this supposed to be the sentinel node? --- NO, sentinel is artificial. In ArborX, it is set to negative 1, where solving between INode and Leaf indexing with the same integer is done by shifting Inode values 'up' the number line, compared to leaf indices, rather than as inverse, like I do with negative values.

2. It appears all INodes either have 2 normal node children or 2 sentinel node children. I am unsure how much of a problem this is.
3. In the INode array, the same child node appears multiple times. In hypotheory, each node should have a single parent, but possibly many 'skip' parents. This seems least like a problem.
    In a single run, only 16 of the 30 leaf nodes appeared as children of the INodes, with 2 duplicates. I could almost be convinced that this is correct, except for the fact that 9 INodes were doubly sentinel. When there should be n-1 INodes for n Leaves. Also, the distribution of INode-to-sentinel should spread across as many INodes as there are levels in the hierarchy. If I got it right, for n Leaves, we expect the following function to find the number of INodes pointing to sentinel:
```julia
function numberoflevels(n)
sum = 0
a = 1
    while a <= n
        a *= 2
        sum += 1
    end
    return println(sum)
end
```
So, that wowuld be 5 levels at 30 leaves. So, we are creating far too many sentinels. Rather, we are not changing them from their original state.

4. In ArborX, they havea convention for deltastar that if i is less than zero or greater than or equal to the the number of internal nodes, then the function returns max values for the numeric type.
    What should my delta do?

    And extending it, for evaluating ranger or rangel with the atomic, I believe I should only ask if rangel/ranger are -1, instead of anything else. The growing convention states that if ranger or rangel are indices marked away from the index of consideration, i, then we reject this evaluation. This change in the code appears to do little, other than remove a single evaluation.

5. Apparent typo from the paper, in line 19 of algorithm 4. They just use delstar of range_left, but in ArborX line 308, it is delstar(range_left - 1). Fixing this issue appears to reduce the average number of untreated INodes from > 10 to ~ 5

6. Their del-star calculates the XOR, and then runs further agumentations based on the index. I am uncertain if this is to reduce divergence by assuming two morton codes are always the same, but it may be more hlepful for me to perform the agumentation than otherwise.
They have:
```
    auto const x = _sorted_morton_codes(i) ^ _sorted_morton_codes(i + 1);

    return x + (!x) * (min_value + (i ^ (i + 1))) - 1;
```
Where '!' is logical NOT, where if x is 0, then the morton codes are the same, and if that value x is zero, it is false, and the opposite of false is true., and i + 1 is my j. Implementing this fix appears to reduce the unchanged nodes even further, now below the anticipated number of sentinels. THe next step is to evaluate the setRope function of ArborX.

## 24 Nov
Still stuck on a few problems. The first is that the root node, supposedly I[1], has a sentinel skip rope. Second, I[n-1] is a double skip rope. Third, some nodes, I and L, are repeated in the ordering, which I don't think is supposed to be possible. Fourth, both L[n] and L[n-1] have sentinel skips, which should only happen to L[n].

According to Figure 1 in Proko, only L[n] should skip to the sentinel. The root should ALSO skip to the sentinel. Also, L[n] should not be discoverable by any internal node parent. It can only be found during traversal by skip rope from another leaf. I[n-1] should have a left child and skip to the sentinel.

Perhaps it will be neough to shift the thread launching by 1, but also the sentinel evaluation by 1. Previously, I would only adjust the thread launch, not the sentinel condition as well!

Making this modification causes multiple solutions to an index of -7, which would be an internal node of index 7, with only 7 leaves. This won't work, but how can it be modified?



I thought I found a mistake, they use a setrope function to and they always +1 the value before application to the Leaf or INode .skip, but in my code it is already set to + 1 in the definition of r. However, I think I am closer, because for the first time, I successfully got a reference to Leaf n without examining all n optoins, except this was done by INode, not by a leaf. There is a plus 1 somewhere that I do not have. WEHRE

I don't know. I have sorted through Proko and ArborX for many hours. I have no idea what the problem is supposed to be at this point. I just do not understand. My structure is stopping 1 unit short, but if I extend it to allow itself to resolve the last unit, it then breaks! Prok's algorithm, as presented in the publication and ArborX should not be capable of creating the tree observed in Figure 1 of their publication. Which perhaps is fair, as it is only an example stackless tree. But still!

## 28 Nov
Okay so I just committed a version in which the as-close-as-possible to the GitHub implementation generates the same result as my version, which took from their publication. This means that my negative indexing to differentiate internal and leaf node does not have a subtle influence on the product, and that I must look elsewhere to explain my problem.

Separately, I realized why nodes were self referencing, and this was a side effect of copying code (from the untested publication pseudo-code) without understanding it. My code experienced 3 bouts of self referencing out of a possible 8 because I incorrectly attatched the left child of the q'th INode as the index number i, instead of the apetrei parent, p. When the entire if else statemtent of ``` del < dell ``` was about evaluating the parentage of the object we are currently considering, i.

So we are down to 3 flavors of issue that may be solved in any number of ways.
1. Not enough INodes are being directed to the sentinel.
2. In fixed position, 3 Inodes skip to the same leaf. In repeated random runs, there does not seem to be a fixed pattern for INode or Leaf repetition, but both infact appear.
3. General correctness issues.

- reversing morton code sort order does not fix (and does some other weird stuff)
- in my equivalent to ```setRope```, evaluating for sentinel nodes, changing evaluation to number of internal nodes does not fix issue


HEY WAIT we absolute did not fix self referencing, we just shifted the value of the self reference by +1. hah! So i have to keep looking in the area of this indexing.

I am trying to solve it manually, but that's difficult because my method is not fully resolved, and the kokkos method is index shift (and slightly mysterious). I suppose I could solve it manually using Kokkos method as a guide and then search for the nuances there.
Running the code by hand,

How does the iterator run?

### Kokkos openmp
https://kokkos.org/kokkos-core-wiki/building.html
https://github.com/arborx/ArborX/wiki/Build
ArborX supports openMP but idk if they support / how they support other parallelization routines like Kokkos::threads or HPX or whatnot.

## 29 Nov
Okay so I tried running manually running the Kokkos approach with my data set, and the sequential nature of it seemed to be very damaging. I processed the first 6 leaves, and then the parent of leaf six,,, was the root. Which means at least the way my code works is exactly the same as my understanding of the Kokkos material itself. And there are no layers of obfuscation between the behavior of my functions and what I think them to be.  Arriving to process the root  as the last thread will just direct the root the last leaf. Now the question is, when I turn on Threads.@threads, do I get a traceable structure for every leaf? And is there a world in which you do not get a traceable structure for every leaf???????? Horrifying.

The structure,, is different I suppose. . . . Except traversal proceeds from the root to -7. Maybe we lack sufficient morton complexity? Well, that was an optimistic thought, but no. Okay, what about a random data set, as this one was kind of crazy.

Something is definitely broken, as it now always traces the root to the 7th dud node. Now what exactly is the matter?

### specializing delta for leaves and internalsf
I believe specializing a delta function for the leaves and the branches is necessary, because if i is greater than or equal to the number of internal nodes, then the function returns typemax, at least in ArborX. But, except depending on some breakthrough in understanding their sorting and index permuting thing, their number of internal nodes is the same number as the index of the right most leaf node. Thus, to a leaf node, it should receive maximum values at or above input values of L[max]. But internal nodes should receive this typemax return at the number of leaf nodes minus one, (also known as the value of max in L[max] in 0-based indexing.).

Now, is this actually the case? Because in this scheme, an internal index would poll a delta between L[max-1] and L[max] as typemax, rather than whatever value it actually is. This feels to me like a loss of data.

### Differentiating between internal node and leaf node
It  turns out, while I can map to the behavior of ArborX quite easily just be replacing internalIndex by just multiplying everything by negative 1, we are still eft out of the town because this negation hides which nodes are internal or leaf nodes, and likely affects the results down the pipes, creating annoying trees. I can't say this is *the* problem, because I spent a while last night manually drawing out a tree with my morton codes and deltas, and I got the same tree generated by my software. Maybe it is a comprehension thing, I understood the same thing from ArborX that I wrote out in my own code, but I have no idea! Regardless, in a test of a 3 leaf tree, placing the root node's left child returned a value of negative 1. That said to me 'hey, the root is self referencing itself'. However, when I jumped to the place in code that evaluates whether the leftChild is a leaf or an internal node, my evaluation said leaf. But the number itself was already negative, causing a false 'self reference'. Furhtermore, if the value were positive 1 instead of negative, then we would have a fully traversable tree.

If we follow ArborX, or specifically its tea-leaves while avoiding looking too closely, then we should have a single array of gridkeys for both the leaves and the internal nodes, and just use indexes,,, correctly,,, to access internal nodes. Afterall, for traversal, there is already a built in comprehension for leaf vs branch, all branches have left children but no leaves have left children. I should not need additional comprehensions on top of negative or positive values. From there, I should have a better chance of evaluating what my problem is in this codebase, and fixing the issues.

Personally, I am happy to learn that part of the problem was the fogginess of my own window into the numbers being crunched.

### internal and leaf nodes as members of the same array
My first problem so far with this scheme is that the data structures are not the same, the leaves have a bunch of morton codes which do not need to be stored with the Gridkey data type, at all. I believe the grid keys and index to atom should be companion arrays that are co-sorted with the structs. But for now, the data inefficiency will be more than fine.

Another crucial issue is that we will have to pop off the internal nodes or sort a specific index length for the morton codes following a refesh of the data.

Okay so the next question is about mutating a value or returning a different value given an input. 

And the next next question is what the hell happened. I solved a 3 leaf structure just find, but then everything collapses at a 4 leaf. And when i tested the 8 leaf, it felt eerily similar to what I was struggling with beforfe.


Okay so it is working really well now. Out of  10 000 runs on a size 7 tree, we have typically fewwer than 1000 untraversable trees. The issue comes in evaluating the root. Sometimes it will connect to the leaf index, that if branch shifted, would result in a valid tree. There is an ambiguity in the code to work out which may be giving us the issue. If we are evaluating which children belong to the root, then rangel and ranger should be 1 and 7 (more generically, 1 and leaves_count). Both delta values here are max value, so we see an odd situation in which two numbers being equal, Julia will randomly decide to pretend one of them is greater than or less than. Let's try to confirm this behavior with a pocket test.It may not work out too well, if the compiler optimizes the values out. Well, we could just include more evaluation . . . .

The pocket test showed that if a and b are the same integer, then ```a < b``` is always false. This false evaluation is mirrored by the same evaluation in the definition of q in which even if we get a bad tree, ```a < b``` is also false. So we have to peruse further.

The issue is firmly an execution order problem. When I start multithreading evaluation from the last leaf towards the first, then we usually acquire the correct tree. But when I run threading in the same order, then usually we aquire the wrong tree. How can this order evaluation issue be ameliorated?

It is an interesting experiment to see if increasing the number of threads solves the issue, but it is possible that increasing the number of threads from 4 to 8 on such a small and brief workload, on a 4 core 8 thread system, does not relieve the evaluation order problem we have, as the threads (leaves) 5 through 7 will be physically evaluated after the first 4.


Okay so, infuriatingly, my modification of leftChild if it was a leaf had a typo in it.So that meant it was a new local variable that went unused. Fixing that issue make the results look like we have regressed all the way to a prior point, in which the root is once again a child, branches hold themselves as their own children, and the tree is nontraversable. In other words, I am reasonably certain that I have simply rewritten the same broken algorithm I had last night, but fixed up to and atmost 1 issue about the print out scores at the end.

## 1 Dec
### traversal batch testing
The results from these batch tests are so strange. When we follow the ArborX implementation by evaluating if branch node inside of the isRightChild code only, we randomly achieve some successful runs. And I have seen the batch swing from 94 successes to 49 successes with the only difference being when I clicked the Run button. And then I decide to remove that logic and toss in a branch_index() call, and that seems to break things a different way by placing large index values in the tree. 

## 6 Dec
### thoughts on threading model
I need to work on my comprehension ofwhat I want the algorithm to achieve before continuing. I also believe that we have to investigate the threading model at hand, as successful tree generation is seemingly enhanced by reverse execution. For certain, there is a runtime / scheduling component to why this does not work. Confoundingly, a batch of 10 000 ~7 leaf tree-generating runs will return dramatically different error and success rates from batch of runs to batch of runs. But it does appear to be the case that reverse order execution favors correct tree formation. It is important to note that either my personal understanding of the algorithm is incorrect, or serial execution will prevent a tree from correctly forming. Curiously, if I iterate for the batch execute, I get some different results that are not inherently clustered, I would say. Well, if this is done by a dirty for-loop at global scope, then yes each result is different. Running 1000 runs of from the position vector to tree traversal 100 total times produces seemingly consistent results: 
If my approach follows arborX and only shifts the value of leftChild inside the isLeftChild=false block, and iterating from L[n] down towards L[1], then about 80% good trees and 20% bad trees.

If I only reverse the thread exuection order, so now it is L[1] iterating to L[n], then about 96% selfish trees and 4% bad trees.

Exchanging Threads.@threads with Threads.@spawn and using backward iteration to 99.5% bad and 0.5% good. Using forward iteration, we get 99.9% bad and 0.1% good. 

Threads.@spawn is apparently nondeterministic in scheduling, so it would be curious to explore structure deformation here.


Changing the scheduling options, in forward iteration, `:dynamic` and `:greedy` are normal, while `:static` changes to a 22.7% bad adn 77.3% selfish. While backwards iteration, `static` changes to 11.4% good and 88.6% bad tree generation. What insanity!

### perf testing on force_coulomb!
Julia's SIMD macro does not do much in any configuration, but I am able to cut execution time on 80 atoms from 4.1 to 3.1 milliseconds by removing bounds checking, and further down to 1.5 ms with the threads macro. Just macros everywhere! Applying these forward to lennardjones and update_pairlist! seemed to do some good.

## 8 Dec
### bounding volume expansion
In ArborX, they generate a local-to-each-thread variable called bounding_volume, that acquires values from the leaf of the current iteration, i.e. from the original value of i. However, in deciding the bounding_volume of a parent, half of the leaves are culled consistent with the thread chopping. I suppose this points to a preprocessing of the bounding volumes or the madness of the skip rope method somehow makes sure it all works in the end. The only unfortunate detail is they have an entire C file devoted to the word and function `expand`. This may take some time, hah!. 

Trying to implement it on my own, I can see a few issues. First struggle is trying to have a flow through structure in which we don't have to evaluate if the current volume in bounding_volume is already inclusive of the expander term. Secondly, I don't want to evaluate if whether the min, max, both, (or maybe neither) of the expander term would actually increase the volume of the bounding)volume variable. Thirdly, I believe it is possible for a parent to only be partially inclusive of all it's children, i.e. that one or more dimensions of any single child may be clipped off due to a boundary child being overall further away in the sum of its x y z dimensions than our single child is any single dimensions. In other words, the bounding volume of the root is not necessarily the entire scene fully encapsulating the centroids of the atoms. Overall, I want to avoid as many instruction and data conditional evaluations as possible and I just want it to flow through.

For starters, the left child of a node is down to the min dimension and the right child is up to the max dimension.


### oh that's smarts
So I believe I finally fixed by broken tree generation. I was setting up for the bounding volume stuff above, and in one crucial place my code differed from ArborX in the labeling of variables, the last if-statement in the while block's ```isLeftChild==false``` block. I originally had:
```julia
if leftChild == rangel
    leftChild = branch_index(leftChild, spec)
end
```


which was then updated to:
```julia
leftChildIsLeaf = (leftChild == rangel)
if !leftChildIsLeaf # HOLY COW THIS WAS WRONG
    leftChild = branch_index(leftChild, spec)
end
```
Further working it, we can remove the second evaluation altogether (though, I am optimistic the compiler would figure this out) by bundling it into the evaluation for bounding volumes.
Fixing that issue has raised the success rate on my tree generation to 100%. So either it is now finally working OR my is_traversable function has logic deficiencies. Either or. Nonetheless, I feel much better. This has been an excruciating lesson in code comprehension

### I want this as a feature drop

An animated traversal for a given point, with appearing and disappearing rectangular prisms as we descend. That would be so cool!!! I have absolutely no idea how to do it haha! But if I get the code working and somehow import a SYCL/LevelZero/CM ET CETERA implementation that can pull data from there into Julia, then we will have a truly interesting reuse situation. Plus, I believe it would be truly important for uhm, verifying the success of my bounding volumes routine.