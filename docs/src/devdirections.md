# Development Directions
   
## Absurd and Obscene Feature Requests
- [] BVH traversal animation as described in Dev Diary, '8 Dec'
- [] Modified `simulate()` run that auto selects algorithms based on user-selectable performance or precision. Algorithms ideally fitting data, but not in a machine learning kind of way. Maybe only applicable to neighbor list algorithms. . .
- [] Multithreading with distributed data that does not have to be reallocated except for chaning data circumstances. Heck, this bullet point should instead be "Determine how redundant allocations are made in opening threads for a single function call and then closing them only to reopen them 1 or 2 function calls later, or if the compiler optimizes this seemingly silly behavior away".

## Roadmap
### Version 0.00.1 - the very beginning

- [x] System instantiation for a generic method and object
- [x] Generic methods, structures, and types to establish an architecture --grtting there
- [x] Separated Testing module
- [x] Generation of a simulation system from generic generator functions
- [x] Generic simulation
- [x] Profiling of a generic simulation and the MDInput?
- [x] Invert generate_object and collect_objects to reduce allocations and prevent having to write messy ass logic that collates elements of object into vectors of object_collection
- [x] Lewrn how to initialize objectCollections in 1 step where they are initialized  and the method of filling is dependent on the type of object, not the collection.
- [x] I think maybe i should not break out the logic for generateposition from collect_objects because doing so might be annoying?
- [x] Threads.@spawn velocityfill position fill ---- that's excessive, no
- [x] Get rid of GenericCollector, toss it into the GenericSystem? No, it's just an initialization. But can pull the initial number of atoms from MDInput


### Version 0.00.2 - foundations

- [x] From initialization of a genericObjectCollection, or the function call of a simulation, generate a pairlist of every unique pair. Hoepfully design this function to be extensible to any number of unique groupings, rather than only pairs -- yes but only with a library
- [x] Vectorize unique pairs, just for shits, giggles,  and curiosity testing? --- no, julian for loops may be faster than vectorized process
- [x] Refactor of GenericObjectCollection and functions to be a single array filled with typed arrays as fields, rather than an object with several attached arrays, so more continuity in data --- no, StructArrays.jl doesnt work or doesnt make sense to me
- [x] Following GOC Refactor,  consider if math should be done on the arrays of GenericObjectCollection, or if the values at start should be copied to smaller arrays that contain only the details pertinent to a given calculation. For now, i think no. It shouldnt be too much of a fuss to create a memory minimal version later, one that copies and cuts out unrelevant information for the sake of parallelization
- [x] Use DataFrames or another package (https://discourse.julialang.org/t/matrix-column-row-labelling/84064/3) for the purpose of assigning labels to the arrays of ObJectCollection, so that code can either call the right subfield by name or by index, rather than index alone and depending on the user to intuit the right name
- [x] Get the docs working
- [x] Get a commentaries/notes doc going in the docs pages. Maybe blog style
- [x] Microbenchmarking of the old genericobjectvollection vs the dataframe one and whether simulate! needs to have the coords and vectors extracted out of it
- [x] Implementstion of a lennard jones potential for generi. Particles
- [x] In simulate!, before the step iterator, we need a section to precalculate values and initialize vectors ---- ??, they are already initialized, weirdo
- [x] To gnericobjectcolelction add radius vector, potential, force
- [x] NaivePairlist algorithm
- [x] Modify the nested logger function by passing it each local variable that it needs to use
- [x] Add a pruner to positions for ones too close at initialization
- [x] Make the pruner time stable
- [x] Set makie work craft into its own module so we dont ask the guthub action to precompile NaiveDynamics when it doesnt have a visualization routine
- [x] Why does ci.yml exist? What does it aim to do? ~run the test package, boss
- [x] Change ci.yml to avoid indicating OS interoperability
- [x] Cutdown on the slop in Simulator
- [x] develop independent methods for MVec  until and unless nonindexable SVec's start winning
- [] improve naming for Vec3D and Stat(ic)?Vec3D
- [x] NaiveLennardJones based on MVec
- [x] NaiveCoulomb based on MVec
- [x] Naive Logging and storage of data as a text file by snapshotting the whole struct~~
- [x] change structuring, so that the object collection is not nested within misc structs
      get rid of the name pile at the start of simulate!()
      simulate takes 3 arguments, a collector, a collection, and the simSpecs, where type dispatch is based around the type of the SimSpecs, if it is a VerletVel or otherwise integrator
      add documentation to describe the arguments for simulate!(), as objectcollection will be shortened to sys, collector to clct and simulationspecification will just be spec
- [x] update simulate!() and object collection so that the force is the force of the current step, so that for a logg of position and force, the listed force sum is the force that (along with velocity and other methods) that caused the particles to change positions between the previous step and the current step
- [x] get rid of dumloop_product!() as it is just an unnecessary composite of larger pieces
- [x] GLMakie integration and MP4 deliverable for data analysis
- [x] add temperature rescaling to catch molecules that suddenly acquire an obscene velocity. Side effect: now only these super fast pairs have any velocity left over
- [x] fix bugs that cause particles to exit the box
- [x] fix bugs that cause unnecessary allocations in VelVer
- [x] notice that a simulation records n_steps + 1 position sets, when trying using frameintervals of 10s, have to shift the value by 1



### Version 0.00.3 - towards a half formal repository
- [] make all Big functions part of a a public API so that they can be tested and developed easier. Especially assembly functions with host functional functions.
- [] fix parametric types in neighborsearch, because the T and K switch positions, replace with I and F? or at least make them consistent
- [] update code naming to reflect the fact that AABB's are only first generated right immediately before bvh traversal, and squash down redundant data structures if at all possible
- [x] fix broken performance by tuple allocation hell, consider switching pairslist to an MVector for values overwrite or trying named tuple shenanigans?
- [] fix velocity rescaling / substitute with alternative method. fix behavior of interactions and parameterization in order to prevent crazy molecular behavior
- [] force LJ may not work correctly. I might have just broken it, but i am uncertain that it correctly calculates the component forces, isntead of just assigning the overall force to each dimension, or something else entirely! ----- TENTATIVELY FIXED, pairslist was messed up generating NaNs and also not doing anything


- [x] fix broken update_pairslist


- [] Improve design of the Logger to be compatible with makie
- [x] Github work flow for a private uhh workspace
- [x] Github based integrations of the code at start and endpoints
- [x] Figure out how to start getting test coverage and using formal unit testing procedures
- [x] Wrap custome types in functions so that a user can call a function and assign labeled arguments (eg "duration=10"), rather than having nameless and ordered fields
- [x] These wrapper functions may also contain side logic for checking inputs are correct as well as the actual logic to be done on the particular system, as shown in Molly.setup
- [] Output logfile with modification of the set up routine to allow the user to add in a place and a type of output, but defaulting to a generic
- [x] Random generation for each component. Check that this works
- [x] Aqua.jl integration that only tests the local package and not every dependency
- [] Consider putting in architecture to read/write data so we can test coverage with fixed values and compare changes with feature development.
- [] make it so push! log only runs at every selected interval, and also make this match the frame interval for makie by having makie take the simSpec as default framerate
- [] change structure definitions in MDINput to be Vec3D instead of Vector{MVector} etc etc
- [x] add kernel abstractions and AMDGPU and oneAPI and CUDA as formal extensions so that they are only precompiled when the script file to use this package includes 'use cuda'.
- [] change all 2 factor ranges to a 2 length tuple
- [] change vectors of structs to be structs of vectors, and add in relevant infrastructure to enable a resort of say the minboundary to change the order in the exact same way of the other elements of the simulation.


- [] check the naive unique pairs function for correctness. I was kinda just throwing stuff at a wall to see if it worked
- [] fix precision selection so the precision can be selected by user exactly one time and is persistent throughout.
- [] Fix position recording so that the simulation can be logged for a user specified number of runs
- [] add ```simulate!()``` resolution so that the system can log the last few steps, if the last step does not trigger a logging of the chunk
- [x] fix bug in simulation recorder where the chunk_index has to be updated inside the for each step loop. when placed inside the record_simulation if statement, then the value will be reset by simulate!() to it's initial definition value each step, even if the place where the value was defined as '2' sits outside of the stepper loop. this could be automagically fixed when we move to more direct function arguments rather than the equivalency pile up top.
- [] fix velocity verlet to prevent velocity from depreciating for no reason. most likely, the velocity values are being overwritten by intermediates, which are based on forces. as forces tend to zero, so shall intermediates and velocities. or the force is just whacked up. not sure!
- [x] use for each map!() for all instances of IntermediateVector = DataVector
- [x] allow record_video() to have user input for the frame recording interval. do this by pushing every multiple of frameInterval to the positions vector
- [x] make these userfill parameters easy to fill in, for name awareness of each paremter
      by having a function of the same name fxn(; param, param, param, defaultparam=1)
information based on other things the user input, like if single precision, then morton encode to 32 bits.
- [x] currently, MakieExt redefines the record_video function stored in PkgExtensions. Will the extension continue to work if it exports record_video on its own?
- [] package extensions methods break upon trying to use them at all because something something Julia doesnt work. In my Dev environment, Iwant as little loaded as possible. Thus, the extensions, but I am tetsint in my dev environment, which means I don't get to use the extension functionality. I believe it would work better for a user situation, in which the Julia environment is not this package's source code. idk
- [] struct instantiate with function for neighborsearch items, so changes to the API are more clear to impelment (but also slightly more tedious)
- [] companion arrays of morton codes, indices to atoms, and grid aabbs and simplified structures for more purposeful datamanagement. These optimizations won't especially work until we have struct arrays and or the deep compression used in contemporary bvh papers.
- [x] if we use the Julia built in environment instead of our own, could we finally have extensions working correctly, so that we are devved into naive dynamics and using the local dev version wiht a napkin test file, while also being abel to use only the dependencies and extensions we want?
- [] functions don't necessarily have to be in order, a function can call a function that is defined physically below it. Use this concept to make the code prettier and better organized.
- [] investigate if other Julian threading routines produce better results. Polyester and OhMyThreads and Dagger come to mind
- [] api.md only has 2 items on it. Why?
- [] fix upper functions of Proko to only iterate over  specific indices of the grid keys array
- [] for boundary expansion, is the for loop flowthrough evaluation effective, or is copyto! effective enough
- [] is it a problem that the root doesnt get updated to cover the whole entire range? may or may not just be a side effect of the algo. ArborX predefines the root.
- [] update ci.yml for a different OS test and to resolve warnings related to chagnes to GHActions


### Version 0.00.4 - feature extensions
- [] Refine functions , ex: sigma6th and sigma12 should be calculated prior to simulation for each unique radius of objects in our objectcollection --- lord willing the compiler will do this at compile time, but i trust nothing and no one.
- [] Improve the boundary_reflect!() in some way to either reduce frequency of checking (pair list), use an aligned array(s) to broadcast that checks in a single statement rather than 6 if statements, convert wall actions into a potential,something else, or all of the above. At leat make it Naive+ rather than just Naive.
- [] Research how boundary conditions are set so as to avoid assessing the value of every particle to see if it exists in the box or not at each time step
- [] optimize naive implementations so that they dont endlessly allocate temporary values and obtain pre-allocated overwrite spaces prior to entering the for-each-step loop
- [] maybe even create a few fancy implementations or Naive+ Naive++ 
- [] integrate tree based neighbor finding
- [] modify the makie extension with the advice posted on their documentation https://docs.makie.org/stable/explanations/recipes#Full-recipes-with-the-@recipe-macro
- [x] fix makie extension so that I don't have to load in GLMakie and all of it's dependencies every single time.
- [] fix up collection.current step and how it is updated inside of simulate!(). it is silly to have to allocate a vector filled with the same data point for each particle at every step. But also, is it really a big deal?
- [] organize helpGwen.md
- [] test out and redevelop struct of arrays of arrays for the Log of ObjectCollections and get a write up on how it's going. it went poorly last time and I am not certain why and I would have to manually search the diary to see if I wrote anything. and maybe i wrote nothign
- [] integrate the julian testing packages as part of a refactor to make naming consistent but also make it easier. for instance, I keep mispelling simLog as simlog when simlog works absolutely fine and syslog might make more sense. or just log.
- [] fix documentation syntax so that documenter.jl transforms the markdown correctly
- [] extend zero() so that it works correctly for a Vec3D and we simpl.ify the zeroing of forces before new calculations
- [] think about how pairwise forces codes have common-boilerplate, can this be abstracted away?
- [] spread package into more files, reduce code weight on each file (.5 instead?)
- [] modify LJ potential of random run to have a sigma and epsilon for each particle
- [] user requestable plots with generic generation method - so we can track the position of particle i throughout a sim, or the mean velocity for a specific duration range and assume these are pictures generated in the local directory as the file that created them
- [] Logging of velocity and position (and any other dynamic property) at a selectable interval
- [] in functions, select whether a CPU, multiCPU, or GPU is to be the analytical device. see how they do this in Molly, as using keyword arguments in function defintions and having different default values cannot be selected for. multiple dispatch only works on types, not on specific fxn inputs.
- [] abstract away force computations to include user defined force weights but also user defined forces.
- [] along with above, user specified interations with type Unions that expect either false, or a parameter. e.g. velocity dampening on a simple rescale is a false on vrescale, or a parameter in the selected Float
      use 'pruning' functions of the type informatino users fill out to make the types consistent, so no multiple Float32(input), figure that all out in the package
- [] test by ony specified tests, rather than the whole package each and every time ichagne a letter or two

#### Style guide things
- [] are all mutation functions inidcated properly?
- [] in mutation functions, is the first operand always the one being mutated?


### Version 0.00.5

- [] Momentum calculations for particles of a selectable and variable radius radius so they bounce against each other for Newtonian-based simulation
- [] Makie rendering / Refactor Makie extension to depict the variable radii of the particles
- [] Measure energy conservation, explore how it evolves
- [] Analyze how to improve the oneloop simulation, write-up in devdiary, investigate why allocation crazy
- [] By lazyarrays?
- [] By ArrayFire, non julian kernel abstraction library?
- [] By replacing GenericObjColl with a vector of Tuples that contain alll of the information? Maintain broadcasting functionality by a vector of tuples of numbers, M/SVectors of numbers, and strings
Version
- [x] Should the component forces LJ, Coulomb etc. be dumped at the very end of each step, given that they are completely recomputed in the next step based on the old situation, rather than additive? -- dumped prior to recalculation
- [] why do particles tend to have velocity almost mostly in the z-direction?
- [] parameterization processing / case study report for relative box size, particle interaction radii and magnitude, temperature, temporal resolution ?
- [] allow user specification of what properties are to be logged
- [] improve update_mortoncode2! to dump all of the bits from a grid value into the morton code at once, if it is possible? there should be some sort of way without iteration to point to every third value, and then just 'or' dump it into inbit, and then just 'or' dump these right into the morton code. Should be abble to get rid of the m loop and only loop for each dimension.
- [] consider working the morton code to go directly from integer coordsinates to lexico sorting, as sucggested on the4 Z-order curve wiki
- [] have a documentation structure that describes our components as Generic, rather than naming all of them Generic. Rename everything with shortened terminology, and find an automated tool to perform this for us? At least, have a Test Everthing button to throw in before and after the edits are made
- [] in line with above, change type dispatch in all of the function to variables rather than fixed structured, if this is wise. I think it is a Normal or Julian thing to do, so to avoid these problems. but maybe not!


### Version 0.01

- [x] Changename from NaiveMD to NaiveDynamics
- [] System initialization from an input file, from a hand constructed input
- [] Definition of a simple particle
- [x] Velocity verlet-based calcuation of stepwise forces, velocities, and positions
- [x] Modeling of spheroids with a lennard jones potential
- [] Particle in a well simulation where the box does something based on the particles being equal or less than a constant distance too close to the wall
- [] Render spheres bouncing against each other in a prism

### Version 0.02

- [] Naive construction of required and assumed unit definitions or importation of unitful.jl for Atom and AtomCollection
- []





## Naive implementations
1. NearestNeighbors  -- calculate distance of every unique and evaluate which are within a threshold distance.

## Fancy implementations