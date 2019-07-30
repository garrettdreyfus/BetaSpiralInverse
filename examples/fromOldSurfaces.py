fileObject = open("data/286364.pickle",'rb')  
surfaces = pickle.load(fileObject)

fileObject = open("data/1500NoNorwegian.pickle",'rb')  
offsets,profiles,deepestindex = pickle.load(fileObject)

surfaces = nstools.convertOldSurfaces(surfaces)

surfaces = nstools.addDataToSurfaces(profiles,surfaces,2)
surfaces = nstools.addStreamFunc(surfaces,profiles)

surfaces =nstools.addXYToSurfaces(surfaces)

interpolatedsurfaces = {}
neighbors={}
lookups={}
for k in surfaces.keys():
    surfaces[k] = nstools.removeDiscontinuities(surfaces[k],radius=0.1)
    interpolatedsurfaces[k],neighbors[k] = nstools.interpolateSurface(surfaces[k])
    lookups[k] = nstools.trueDistanceLookup(interpolatedsurfaces[k],neighbors[k])

graph.graphSurfaces(surfaces,"s")
staggeredsurfaces = nstools.addHorizontalGrad(interpolatedsurfaces,neighbors,lookups)
staggeredsurfaces = nstools.addHAndVerticalGrad(staggeredsurfaces)
staggeredsurfaces = nstools.addK(staggeredsurfaces,"data/bathVar.pickle")

