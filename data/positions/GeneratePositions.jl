using JLD2
using NaiveDynamics

thresh = 0.03

clct10 = GenericRandomCollector(; floattype=Float32,
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
clxn10 = collect_objects(clct10)

clct100 = GenericRandomCollector(; floattype=Float32,
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
clxn100 = collect_objects(clct100)

clct1000 = GenericRandomCollector(; floattype=Float32,
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
clxn1000 = collect_objects(clct1000)

clct2000 = GenericRandomCollector(; floattype=Float32,
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
clxn2000 = collect_objects(clct2000)

clct3000 = GenericRandomCollector(; floattype=Float32,
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
clxn3000 = collect_objects(clct3000)

clct5000 = GenericRandomCollector(; floattype=Float32,
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
clxn5000 = collect_objects(clct5000)

clct20000 = GenericRandomCollector(; floattype=Float32,
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
clxn20000 = collect_objects(clct20000)

pos10 = clxn10.position
pos100 = clxn100.position
pos1000 = clxn1000.position
pos2000 = clxn2000.position
pos3000 = clxn3000.position
pos5000 = clxn5000.position
pos20000 = clxn20000.position

jldsave("data/positions/positions.jld2"; pos10, pos100, pos1000, pos2000, pos3000, pos5000, pos20000 )
