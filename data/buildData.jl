module buildData
using NaiveDynamics


# 1. Build a large simulation with set values.
valuesCollector = GenericUserValueCollector(100, 0.2, 12.9)
valuesCollection = collect_objects(valuesCollector)
valuesSystem = GenericSystem(10, 1, 1, valuesCollection)
valuesSimulation = GenericSimulation(valuesSystem, true)
simulate!(valuesSimulation)

# 2. TODO build a large simulation with random values (will take a while to run)

end