##  not sure yet what i want to test here.
## that position spawning will never runaway? eh.
## i mayt want to force interesting spawning patterns with different distances tricks
# myCollector = GenericRandomCollector(; floattype=Float32,
#                                     objectnumber=8,
#                                     minDim=tuple(0.0, 0.0, 0.0),
#                                     maxDim=tuple(1.0, 1.0, 1.0),
#                                     temperature=0.01,
#                                     randomvelocity=false,
#                                     minmass=1.0,
#                                     maxmass=5.0,
#                                     minimumdistance=1.0,
#                                     mincharge=-1f-9,
#                                     maxcharge=1f-9
# )

# @testset "position spawning" begin
#     @test myCollection = collect_objects(myCollector)
# end