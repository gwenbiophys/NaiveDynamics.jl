using JLD2
using NaiveDynamics

thresh = 0.03

clct10 = GenericStaticRandomCollector(; floattype=Float32,
                                objectnumber=10,
                                minDim=tuple(0.0, 0.0, 0.0),
                                maxDim=tuple(1.0, 1.0, 1.0),
                                temperature=0.01,
                                randomvelocity=false,
                                minmass=1.0,
                                maxmass=5.0,
                                minimumdistance=0.0001,
                                mincharge=-1f-9,
                                maxcharge=1f-9
)

clct100 = GenericStaticRandomCollector(; floattype=Float32,
                                objectnumber=100,
                                minDim=tuple(0.0, 0.0, 0.0),
                                maxDim=tuple(1.0, 1.0, 1.0),
                                temperature=0.01,
                                randomvelocity=false,
                                minmass=1.0,
                                maxmass=5.0,
                                minimumdistance=0.0001,
                                mincharge=-1f-9,
                                maxcharge=1f-9
)

clct1000 = GenericStaticRandomCollector(; floattype=Float32,
                                objectnumber=1000,
                                minDim=tuple(0.0, 0.0, 0.0),
                                maxDim=tuple(1.0, 1.0, 1.0),
                                temperature=0.01,
                                randomvelocity=false,
                                minmass=1.0,
                                maxmass=5.0,
                                minimumdistance=0.0001,
                                mincharge=-1f-9,
                                maxcharge=1f-9
)

clct2000 = GenericStaticRandomCollector(; floattype=Float32,
                                objectnumber=2000,
                                minDim=tuple(0.0, 0.0, 0.0),
                                maxDim=tuple(1.0, 1.0, 1.0),
                                temperature=0.01,
                                randomvelocity=false,
                                minmass=1.0,
                                maxmass=5.0,
                                minimumdistance=0.0001,
                                mincharge=-1f-9,
                                maxcharge=1f-9
)

clct3000 = GenericStaticRandomCollector(; floattype=Float32,
                                objectnumber=3000,
                                minDim=tuple(0.0, 0.0, 0.0),
                                maxDim=tuple(1.0, 1.0, 1.0),
                                temperature=0.01,
                                randomvelocity=false,
                                minmass=1.0,
                                maxmass=5.0,
                                minimumdistance=0.0001,
                                mincharge=-1f-9,
                                maxcharge=1f-9
)

clct5000 = GenericStaticRandomCollector(; floattype=Float32,
                                objectnumber=5000,
                                minDim=tuple(0.0, 0.0, 0.0),
                                maxDim=tuple(1.0, 1.0, 1.0),
                                temperature=0.01,
                                randomvelocity=false,
                                minmass=1.0,
                                maxmass=5.0,
                                minimumdistance=0.0001,
                                mincharge=-1f-9,
                                maxcharge=1f-9
)

clct20000 = GenericStaticRandomCollector(; floattype=Float32,
                                objectnumber=20000,
                                minDim=tuple(0.0, 0.0, 0.0),
                                maxDim=tuple(1.0, 1.0, 1.0),
                                temperature=0.01,
                                randomvelocity=false,
                                minmass=1.0,
                                maxmass=5.0,
                                minimumdistance=0.0001,
                                mincharge=-1f-9,
                                maxcharge=1f-9
)

clct100000 = GenericStaticRandomCollector(; floattype=Float32,
                                objectnumber=100000,
                                minDim=tuple(0.0, 0.0, 0.0),
                                maxDim=tuple(1.0, 1.0, 1.0),
                                temperature=0.01,
                                randomvelocity=false,
                                minmass=1.0,
                                maxmass=5.0,
                                minimumdistance=0.0001,
                                mincharge=-1f-9,
                                maxcharge=1f-9
)

clct500000 = GenericStaticRandomCollector(; floattype=Float32,
                                objectnumber=500000,
                                minDim=tuple(0.0, 0.0, 0.0),
                                maxDim=tuple(1.0, 1.0, 1.0),
                                temperature=0.01,
                                randomvelocity=false,
                                minmass=1.0,
                                maxmass=5.0,
                                minimumdistance=0.0001,
                                mincharge=-1f-9,
                                maxcharge=1f-9
)

clct1000000 = GenericStaticRandomCollector(; floattype=Float32,
                                objectnumber=1000000,
                                minDim=tuple(0.0, 0.0, 0.0),
                                maxDim=tuple(1.0, 1.0, 1.0),
                                temperature=0.01,
                                randomvelocity=false,
                                minmass=1.0,
                                maxmass=5.0,
                                minimumdistance=0.0001,
                                mincharge=-1f-9,
                                maxcharge=1f-9
)


pos10 = generate_positions(clct10)
pos100 = generate_positions(clct100)
pos1000 = generate_positions(clct1000)
pos2000 = generate_positions(clct2000)
pos3000 = generate_positions(clct3000)
pos5000 = generate_positions(clct5000)
pos20000 = generate_positions(clct20000)
pos100000 = generate_positions(clct100000)
pos500000 = generate_positions(clct500000)
pos1000000 = generate_positions(clct1000000)

jldsave("assets/positions/positions.jld2"; pos10, pos100, pos1000, pos2000, pos3000, pos5000, pos20000, pos100000, pos500000, pos1000000 )
