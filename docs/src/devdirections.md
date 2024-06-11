# Development Directions

<details>
   <summary>Version wishlist</summary>
    ## Roadmap
    ### Version 0.00.1
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


    ### Version 0.00.2
    - [x] From initialization of a genericObjectCollection, or the function call of a simulation, generate a pairlist of every unique pair. Hoepfully design this function to be extensible to any number of unique groupings, rather than only pairs -- yes but only with a library
    - [x] Vectorize unique pairs, just for shits, giggles,  and curiosity testing? --- no, julian for loops may be faster than vectorized process
    - [x] Refactor of GenericObjectCollection and functions to be a single array filled with typed arrays as fields, rather than an object with several attached arrays, so more continuity in data --- no, StructArrays.jl doesnt work or doesnt make sense to me
    - [x] Following GOC Refactor,  consider if math should be done on the arrays of GenericObjectCollection, or if the values at start should be copied to smaller arrays that contain only the details pertinent to a given calculation. For now, i think no. It shouldnt be too much of a fuss to create a memory minimal version later, one that copies and cuts out unrelevant information for the sake of parallelization
    - [x] Use DataFrames or another package (https://discourse.julialang.org/t/matrix-column-row-labelling/84064/3) for the purpose of assigning labels to the arrays of ObJectCollection, so that code can either call the right subfield by name or by index, rather than index alone and depending on the user to intuit the right name
    - [x] (Get the docs working)
    - [] (Get test coverage working and automated with each commit)
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
    - [] Why does ci.yml exist? What does it aim to do?
    - [x] Change ci.yml to avoid indicating OS interoperability
    - [] Cutdown on the slop in Simulator
    - [] develop independent methods for MVec  until and unless nonindexable SVec's start winning
    - [] improve naming for Vec3D and Stat(ic)?Vec3D
    - [] NaiveLennardJones based on MVec
    - [] NaiveCoulomb based on MVec
    - [] Naive bonding groups
    - [] Naive Logging and CSV block-streaming with MVec 


    ### Version 0.00.3
    - [] GLMakie integration and MP4 deliverable for data analysis
    - [] Improve design of the Logger to be compatible with makie
    - [] Add ElectronDisplays and configure to display plots in their own window
    - [x] Github work flow for a private uhh workspace
    - [x] Github based integrations of the code at start and endpoints
    - [] Figure out how to start getting test coverage and using formal unit testing procedures
    - [] Wrap custome types in functions so that a user can call a function and assign labeled arguments (eg "duration=10"), rather than having nameless and ordered fields
    - [] (These wrapper functions may also contain side logic for checking inputs are correct as well as the actual logic to be done on the particular system, as shown in Molly.setup)
    - [] Output logfile with modification of the set up routine to allow the user to add in a place and a type of output, but defaulting to a generic
    - [x] Random generation for each component. Check that this works
    - [] Aqua.jl
    - [] Consider putting in architecture to read data from input files so we can test coverage with fixed values and analyze for changes with feature development.
    - [] add several Naive implementations
    - [] For instance, sigma6th and sigma12 should be calculated prior to simulation for each unique radius of objects in our objectcollection --- lord willing the compiler will do this at compile time, but i trust nothing and no one.
    - [] check the naive unique pairs function for correctness. I was kinda just throwing stuff at a wall to see if it worked


    ### Version 0.00.4
    - [x] Construction of a solid wall that inverts the velocity term of any particle with a >= value as the wall
    - [] Research how boundary conditions are set so as to avoid assessing the value of every particle to see if it exists in the box or not at each time step
    - [] optimize naive implementations so that they dont endlessly allocate temporary values and obtain pre-allocated overwrite spaces prior to entering the for-each-step loop
    - [] maybe even create a few fancy implementations 
    - [] integrate tree based neighbor finding

    ### Version 0.00.5
    - [] Momentum calculations for particles of a selectable and variable radius radius so they bounce against each other for Newtonian-based simulation
    - [] Makie rendering / Refactor Makie extension to depict the variable radii of the particles
    - [] Measure energy conservation, explore how it evolves
    - [] Analyze how to improve the oneloop simulation, write-up in devdiary, investigate why allocation crazy
    - [] By lazyarrays?
    - [] By ArrayFire, non julian kernel abstraction library?
    - [] By replacing GenericObjColl with a vector of Tuples that contain alll of the information? Maintain broadcasting functionality by a vector of tuples of numbers, M/SVectors of numbers, and strings
    Version

    ### Version 0.01
    - [x] Changename from NaiveMD to NaiveDynamics
    - [] System initialization from an input file, from a hand constructed input
    - [] Definition of a simple particle
    - [x] Velocity verlet-based calcuation of stepwise forces, velocities, and positions
    - [x] Modeling of spheroids with a lennard jones potential
    - [] Logging of velocity and position (and any other dynamic property) at a selectable interval
    - [] Particle in a well simulation where the box does something based on the particles being equal or less than a constant distance too close to the wall
    - [] Render spheres bouncing against each other in a prism

    ### Version 0.02
    - [] Naive construction of required and assumed unit definitions or importation of unitful.jl for Atom and AtomCollection
    - []
</details>

## Naive implementations

## Fancy implementations

## Optimization-type TODO

## Miscellaneous TODO



## Interface of presently implemented functions
```@index
Order   = - [:module, :type, :constant, :function, :macro]
```

```@autodocs
Modules = - [NaiveDynamics]
Order   = - [:module, :type, :constant, :function, :macro]
```