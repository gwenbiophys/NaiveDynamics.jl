module Testing
include("MDInput.jl")
#include("definitionsAndNotes.jl")
include("Simulator.jl")

# 1. Testing GenericObject
"""
dorg= [[0.1, 0.5] [0.2, 0.6] [0.3, 0.7]]
dorg = GenericObject("Duck", [[0.1] [0.2] [0.3]], [[0.4] [0.5] [0.6]], ID() )
dorf = GenericObject("Rabbit", [[0.1] [0.2] [0.3]], [[0.4] [0.5] [0.6]], ID() )
println(dorg, " ", typeof(dorg))
println(dorf, " ", typeof(dorf))
println(voidobject)
"""

# 2. Testing generate_positions!
"""
testCollection = Array{AbstractFloat}(undef, 3, NUMBER_OF_OBJECTS)
println(testCollection)
#sizehint!(testCollection, 10)
generate_positions3!(testCollection)
println(testCollection)
"""

end # module